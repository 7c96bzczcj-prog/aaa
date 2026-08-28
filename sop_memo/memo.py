"""贴在屏幕右上角的置顶便签：TGF-β 撤药模型 SOP（NK-92）第一轮日程。

跑法：
    python sop_memo/memo.py                     # Day 0 = 下一个周二
    python sop_memo/memo.py --day0 2026-09-01   # 指定 Day 0
    python sop_memo/memo.py --decorated         # 万一无边框窗口在你系统上不听话
    python sop_memo/memo.py --test-notify       # 试一条 macOS 通知，确认系统放行了

只用标准库 tkinter。勾选状态与窗口位置存在 ~/.tgfb_sop_memo.json，重开不丢。
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
import tkinter as tk
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import notify  # noqa: E402
from schedule import (  # noqa: E402
    PREP_TASKS,
    ROUND2_GATES,
    TRACK_LABEL,
    TRACK_MAIN,
    WEEKDAY_CN,
    build_round1,
    find_next,
    find_today,
    flow,
    heads_up_line,
    next_weekday,
    notification_for,
    today_line,
    tomorrow_line,
)

STATE_PATH = Path.home() / ".tgfb_sop_memo.json"

# 便签纸配色，浅色底，长时间盯着不刺眼
BG = "#fdf6d8"
BG_HEAD = "#f2e3a3"
BG_TODAY = "#fff2b0"
BG_ROW = "#fdf6d8"
BG_ROW_ALT = "#f8efcd"
BG_AHEAD = "#ffe08a"
FG = "#3a3226"
FG_DIM = "#8a7f68"
FG_DONE = "#a9a08c"
BORDER = "#d8c88a"

MARGIN = 12  # 距屏幕右上角的留白
WIDTH = 400
MAX_HEIGHT_RATIO = 0.88  # 最高不超过屏幕高度的这个比例

TICK_MS = 300_000  # 5 分钟重画并检查一次通知
DEFAULT_NOTIFY_AT = _dt.time(9, 0)  # 每天几点往通知中心发那一条


def pick_font() -> str:
    """挑一个各平台都有中文的字体。"""
    if sys.platform == "darwin":
        return "PingFang SC"
    if sys.platform.startswith("win"):
        return "Microsoft YaHei"
    return "Noto Sans CJK SC"


def load_state() -> dict:
    try:
        with open(STATE_PATH, encoding="utf-8") as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def save_state(state: dict) -> None:
    try:
        tmp = STATE_PATH.with_suffix(".json.tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(state, fh, ensure_ascii=False, indent=1)
        os.replace(tmp, STATE_PATH)
    except OSError:
        pass  # 存不下就算了，不能因此弄崩便签


class Memo:
    def __init__(self, root: tk.Tk, day0: _dt.date, decorated: bool = False,
                 notify_at: _dt.time | None = DEFAULT_NOTIFY_AT) -> None:
        self.root = root
        self.day0 = day0
        self.decorated = decorated
        self.notify_at = notify_at
        self.plans = build_round1(day0)
        self.state = load_state()
        self.state.setdefault("checked", {})
        self.collapsed = bool(self.state.get("collapsed", False))
        self.font = pick_font()
        self._drag = (0, 0)

        root.title("SOP 备忘")
        root.configure(bg=BG)
        root.attributes("-topmost", True)
        if not decorated:
            root.overrideredirect(True)

        self.outer = tk.Frame(root, bg=BG, highlightthickness=1, highlightbackground=BORDER)
        self.outer.pack(fill="both", expand=True)

        self._build_header()
        self.body = tk.Frame(self.outer, bg=BG)
        self.body.pack(fill="both", expand=True)

        self.render()
        self.place_window()
        self.maybe_notify()  # 开机晚了也补发当天那条
        self.root.after(TICK_MS, self._tick)  # 定期重画，跨过午夜自动换“今天”

    # ---------------- 布局 ----------------

    def _build_header(self) -> None:
        head = tk.Frame(self.outer, bg=BG_HEAD)
        head.pack(fill="x")

        left = tk.Frame(head, bg=BG_HEAD)
        left.pack(side="left", fill="x", expand=True, padx=(10, 0), pady=6)
        tk.Label(
            left, text="TGF-β 撤药模型 SOP · NK-92", bg=BG_HEAD, fg=FG,
            font=(self.font, 11, "bold"), anchor="w",
        ).pack(fill="x")
        tk.Label(
            left,
            text=f"第一轮　Day 0 = {self.day0:%Y-%m-%d} 周{WEEKDAY_CN[self.day0.weekday()]}",
            bg=BG_HEAD, fg=FG_DIM, font=(self.font, 9), anchor="w",
        ).pack(fill="x")

        btns = tk.Frame(head, bg=BG_HEAD)
        btns.pack(side="right", padx=6)
        self.toggle_btn = self._flat_button(btns, "▾" if self.collapsed else "▴", self.toggle)
        self.toggle_btn.pack(side="left")
        self._flat_button(btns, "✕", self.root.destroy).pack(side="left")

        # 无边框窗口自己实现拖动
        for widget in (head, left, *left.winfo_children()):
            widget.bind("<Button-1>", self._drag_start)
            widget.bind("<B1-Motion>", self._drag_move)
            widget.bind("<ButtonRelease-1>", self._drag_end)

    def _flat_button(self, parent: tk.Widget, text: str, command) -> tk.Label:
        lbl = tk.Label(parent, text=text, bg=BG_HEAD, fg=FG_DIM,
                       font=(self.font, 11), padx=6, cursor="hand2")
        lbl.bind("<Button-1>", lambda _e: command())
        lbl.bind("<Enter>", lambda _e: lbl.configure(fg=FG))
        lbl.bind("<Leave>", lambda _e: lbl.configure(fg=FG_DIM))
        return lbl

    SCROLL_SEQUENCES = ("<MouseWheel>", "<Button-4>", "<Button-5>")

    def render(self) -> None:
        # bind_all 挂在整个应用上，不随 canvas 一起销毁；收起后再滚滚轮
        # 就会打到已销毁的 canvas 上。重画前先摘干净。
        for seq in self.SCROLL_SEQUENCES:
            self.root.unbind_all(seq)
        for child in self.body.winfo_children():
            child.destroy()
        self.canvas = None
        self.inner = None

        today = _dt.date.today()
        self.focus_widget = None

        self._flow_panel(self.body, today)

        if self.collapsed:
            self._render_focus(self.body, today)
        else:
            self._render_full(self.body, today)

    def _flow_panel(self, parent: tk.Widget, today: _dt.date) -> None:
        """随时间流动的顶部面板：今天做什么、明天做什么，空档长就提前预告。"""
        f = flow(self.plans, today)
        panel = tk.Frame(parent, bg=BG_TODAY)
        panel.pack(fill="x")

        tk.Label(
            panel, text=today_line(f), bg=BG_TODAY,
            fg=FG if f.today_plan else FG_DIM, font=(self.font, 10, "bold"),
            anchor="w", padx=10, wraplength=WIDTH - 24, justify="left",
        ).pack(fill="x", pady=(6, 0))

        tk.Label(
            panel, text=tomorrow_line(f), bg=BG_TODAY,
            fg=FG if f.tomorrow_plan else FG_DIM, font=(self.font, 9),
            anchor="w", padx=10, wraplength=WIDTH - 24, justify="left",
        ).pack(fill="x", pady=(1, 6))

        ahead = heads_up_line(f)
        if ahead:
            tk.Label(
                panel, text=ahead, bg=BG_AHEAD, fg=FG, font=(self.font, 9, "bold"),
                anchor="w", padx=10, pady=4, wraplength=WIDTH - 24, justify="left",
            ).pack(fill="x")

    def _render_focus(self, parent: tk.Widget, today: _dt.date) -> None:
        """收起状态：只显示今天，今天没事就显示下一次。"""
        plan = find_today(self.plans, today) or find_next(self.plans, today)
        if plan is None:
            tk.Label(parent, text="第一轮已结束", bg=BG, fg=FG_DIM,
                     font=(self.font, 10), padx=10, pady=8).pack(fill="x")
            return
        self._day_block(parent, plan, today, alt=False)

    def _render_full(self, parent: tk.Widget, today: _dt.date) -> None:
        canvas = tk.Canvas(parent, bg=BG, highlightthickness=0, width=WIDTH)
        scroll = tk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        inner = tk.Frame(canvas, bg=BG)

        inner.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        window = canvas.create_window((0, 0), window=inner, anchor="nw", width=WIDTH - 2)
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))
        self.canvas, self.inner = canvas, inner

        for seq, delta in ((4, 120), (5, -120)):  # X11 滚轮
            canvas.bind_all(f"<Button-{seq}>", lambda _e, d=delta: canvas.yview_scroll(-d // 120, "units"))
        canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(-1 if e.delta > 0 else 1, "units"))

        prep = self._section(inner, "开始前（SOP 未指定日期，自行安排）")
        # Day 0 之前，该做的就是这些准备，视图先停在这儿
        if self.plans and today < self.plans[0].date:
            self.focus_widget = prep
        for idx, task in enumerate(PREP_TASKS):
            self._task_row(inner, key=f"prep:{idx}", task=task, alt=idx % 2 == 1, dim=True)

        for i, plan in enumerate(self.plans):
            self._day_block(inner, plan, today, alt=i % 2 == 1)

        self._section(inner, "第二轮：窄通道检验（不排日期，条件性启动）")
        for gate in ROUND2_GATES:
            self._note(inner, f"· {gate}")

    def _day_block(self, parent: tk.Widget, plan, today: _dt.date, alt: bool) -> None:
        is_today = plan.date == today
        past = plan.date < today
        head_bg = BG_TODAY if is_today else (BG_ROW_ALT if alt else BG_ROW)
        head = tk.Frame(parent, bg=head_bg)
        head.pack(fill="x", pady=(8, 0))

        mark = "● " if is_today else ""
        tk.Label(
            head, text=f"{mark}{plan.label}　{plan.date_label}",
            bg=head_bg, fg=FG_DIM if past else FG, font=(self.font, 10, "bold"),
            anchor="w", padx=10, pady=3,
        ).pack(fill="x")

        # 视图跟着时间走：滚到今天；今天没事就滚到下一个操作日
        if self.focus_widget is None and plan.date >= today:
            self.focus_widget = head

        for idx, task in enumerate(plan.tasks):
            self._task_row(parent, key=f"{plan.rel_day}:{idx}", task=task,
                           alt=idx % 2 == 1, highlight=is_today, dim=past)

    def _task_row(self, parent: tk.Widget, key: str, task, alt: bool,
                  highlight: bool = False, dim: bool = False) -> None:
        bg = BG_TODAY if highlight else (BG_ROW_ALT if alt else BG_ROW)
        row = tk.Frame(parent, bg=bg)
        row.pack(fill="x")

        done = bool(self.state["checked"].get(key))
        box = tk.Label(row, text="☑" if done else "☐", bg=bg,
                       fg=FG_DONE if done else FG, font=(self.font, 11),
                       padx=8, cursor="hand2")
        box.pack(side="left", anchor="n", pady=1)

        tag = "" if task.track == TRACK_MAIN else f"[{TRACK_LABEL[task.track]}] "
        text = tk.Label(
            row, text=f"{tag}{task.text}", bg=bg,
            fg=FG_DONE if done else (FG_DIM if dim else FG),
            font=(self.font, 9, "overstrike" if done else "normal"),
            anchor="w", justify="left", wraplength=WIDTH - 96, cursor="hand2",
        )
        text.pack(side="left", fill="x", expand=True, pady=1)
        tk.Label(row, text=task.ref, bg=bg, fg=FG_DIM,
                 font=(self.font, 8), padx=6).pack(side="right", anchor="n", pady=2)

        for widget in (box, text):
            widget.bind("<Button-1>", lambda _e, k=key: self.toggle_task(k))

    def _section(self, parent: tk.Widget, title: str) -> tk.Label:
        lbl = tk.Label(parent, text=title, bg=BG, fg=FG_DIM, font=(self.font, 9, "bold"),
                       anchor="w", padx=10)
        lbl.pack(fill="x", pady=(10, 2))
        return lbl

    def _note(self, parent: tk.Widget, text: str) -> None:
        tk.Label(parent, text=text, bg=BG, fg=FG_DIM,
                 font=(self.font, 9), anchor="w", justify="left",
                 wraplength=WIDTH - 28, padx=10, pady=3).pack(fill="x")

    # ---------------- 交互 ----------------

    def toggle_task(self, key: str) -> None:
        checked = self.state["checked"]
        if checked.get(key):
            checked.pop(key, None)
        else:
            checked[key] = True
        save_state(self.state)
        self.render()

    def toggle(self) -> None:
        self.collapsed = not self.collapsed
        self.state["collapsed"] = self.collapsed
        save_state(self.state)
        self.toggle_btn.configure(text="▾" if self.collapsed else "▴")
        self.render()
        self.place_window()

    def place_window(self) -> None:
        """贴到屏幕右上角；若之前拖动过，就用记住的位置。"""
        self.root.update_idletasks()
        screen_w = self.root.winfo_screenwidth()
        screen_h = self.root.winfo_screenheight()
        budget = int(screen_h * MAX_HEIGHT_RATIO)

        # Canvas 自己不申请高度，得按内容算一次，超出屏幕预算的部分交给滚动条
        if self.canvas is not None and self.inner is not None:
            chrome = self.root.winfo_reqheight() - self.canvas.winfo_reqheight()
            self.canvas.configure(height=max(160, min(self.inner.winfo_reqheight(), budget - chrome)))
            self.root.update_idletasks()

        self._scroll_to_focus()
        height = min(self.root.winfo_reqheight(), budget)

        pos = self.state.get("pos")
        if isinstance(pos, list) and len(pos) == 2:
            x, y = int(pos[0]), int(pos[1])
            x = max(0, min(x, screen_w - WIDTH))
            y = max(0, min(y, screen_h - 80))
        else:
            x, y = screen_w - WIDTH - MARGIN, MARGIN
        self.root.geometry(f"{WIDTH}x{height}+{x}+{y}")

    def _scroll_to_focus(self) -> None:
        """把今天（或下一个操作日）滚到可视区顶部，过去的留在上面不挡路。"""
        if self.canvas is None or self.inner is None or self.focus_widget is None:
            return
        self.root.update_idletasks()
        total = self.inner.winfo_height()
        if total <= 0:
            return
        self.canvas.yview_moveto(max(0.0, (self.focus_widget.winfo_y() - 6) / total))

    def _drag_start(self, event: tk.Event) -> None:
        self._drag = (event.x_root - self.root.winfo_x(), event.y_root - self.root.winfo_y())

    def _drag_move(self, event: tk.Event) -> None:
        x = event.x_root - self._drag[0]
        y = event.y_root - self._drag[1]
        self.root.geometry(f"+{x}+{y}")

    def _drag_end(self, _event: tk.Event) -> None:
        self.state["pos"] = [self.root.winfo_x(), self.root.winfo_y()]
        save_state(self.state)

    # ---------------- 通知中心 ----------------

    def maybe_notify(self, now: _dt.datetime | None = None) -> bool:
        """到点就往通知中心发一条。每天最多一条，发过就记下日期。

        没设 ``--notify-at``、还没到点、当天已发过、或今天没什么可说的，都不发。
        """
        if self.notify_at is None:
            return False
        now = now or _dt.datetime.now()
        today = now.date()
        if self.state.get("last_notified") == today.isoformat():
            return False
        if now.time() < self.notify_at:
            return False
        note = notification_for(flow(self.plans, today))
        if note is None:
            # 今天没内容也记一笔，免得每 5 分钟重算一次
            self.state["last_notified"] = today.isoformat()
            save_state(self.state)
            return False
        sent = notify.send(note.title, note.subtitle, note.body)
        self.state["last_notified"] = today.isoformat()
        save_state(self.state)
        return sent

    def _tick(self) -> None:
        self.render()
        self.place_window()
        self.maybe_notify()
        self.root.after(TICK_MS, self._tick)


def parse_time(text: str) -> _dt.time:
    """把 ``HH:MM`` 解析成时间，给 argparse 用。"""
    try:
        return _dt.datetime.strptime(text, "%H:%M").time()
    except ValueError:
        raise argparse.ArgumentTypeError(f"时间要写成 HH:MM，收到的是 {text!r}") from None


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="右上角 SOP 备忘便签")
    parser.add_argument("--day0", help="Day 0 日期 YYYY-MM-DD；缺省为下一个周二")
    parser.add_argument("--decorated", action="store_true", help="保留系统窗口边框")
    parser.add_argument("--reset", action="store_true", help="清空勾选与记住的位置")
    parser.add_argument("--notify-at", type=parse_time, default=DEFAULT_NOTIFY_AT,
                        metavar="HH:MM", help="每天几点发通知中心提醒，缺省 09:00")
    parser.add_argument("--no-notify", action="store_true", help="完全不发通知")
    parser.add_argument("--test-notify", action="store_true",
                        help="立刻发一条测试通知然后退出，用来确认系统放行了通知")
    return parser.parse_args(argv)


def run_test_notify() -> int:
    """发一条测试通知。macOS 上通知权限是最常见的失败点，给个当场能试的法子。"""
    if not notify.available():
        print(f"当前系统（{sys.platform}）不支持通知中心，便签本身照常用。")
        return 1
    ok = notify.send("TGF-β SOP · 测试", "通知中心已接通", "看到这条就说明提醒能发出来。")
    if ok:
        print("已发出。没看见的话到 系统设置 → 通知 → 脚本编辑器 里把通知打开。")
        return 0
    print("发送失败。检查 osascript 是否可用，以及系统设置里的通知权限。")
    return 1


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.test_notify:
        return run_test_notify()
    if args.reset:
        STATE_PATH.unlink(missing_ok=True)
    day0 = _dt.date.fromisoformat(args.day0) if args.day0 else next_weekday(_dt.date.today(), 1)
    root = tk.Tk()
    Memo(root, day0, decorated=args.decorated,
         notify_at=None if args.no_notify else args.notify_at)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
