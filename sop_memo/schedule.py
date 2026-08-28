"""TGF-β 撤药模型 SOP（NK-92）v1.2 的日程计算。

全部内容照抄 SOP 的相对天安排（Day 0 = 铺板当天），只把相对天换算成日历日期。
纯逻辑，不依赖 tkinter，便于单独测试。
"""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field

SOP_TITLE = "TGF-β 撤药模型 SOP（NK-92）v1.2"

WEEKDAY_CN = ["一", "二", "三", "四", "五", "六", "日"]

# 主实验轨 / 并行耐受板轨 / 每个取样日的固定记录动作
TRACK_MAIN = "main"
TRACK_TOLERANCE = "tolerance"
TRACK_RECORD = "record"

TRACK_LABEL = {
    TRACK_MAIN: "主板",
    TRACK_TOLERANCE: "耐受板",
    TRACK_RECORD: "记录",
}


@dataclass(frozen=True)
class Task:
    text: str
    ref: str
    track: str = TRACK_MAIN


@dataclass
class DayPlan:
    rel_day: int
    date: _dt.date
    tasks: list[Task] = field(default_factory=list)

    @property
    def label(self) -> str:
        return f"Day {self.rel_day}"

    @property
    def date_label(self) -> str:
        return f"{self.date:%m-%d} 周{WEEKDAY_CN[self.date.weekday()]}"


# ---------------------------------------------------------------------------
# 第一轮：撤药模型建立（SOP §5）
# ---------------------------------------------------------------------------

# SOP §5.2 时间表 + §5.5 取样 + §5.6 并行耐受板 + §9 记录
_ROUND1: dict[int, list[Task]] = {
    0: [
        Task("计数并调至 3×10⁵/mL，1 mL/孔铺 24 孔板", "§5.1", TRACK_MAIN),
        Task("按组加药：对照 IL-2；持续/撤药组 +TGF-β 10 ng/mL", "§4.2 §5.1", TRACK_MAIN),
        Task("非 TGF-β 组同步加等体积酸性 BSA 载体（4 mM HCl + 0.1% BSA）", "§4.2", TRACK_MAIN),
        Task("取 Day 0 样本 → 立即流式染色", "§5.5 §6", TRACK_MAIN),
        Task("另开耐受板：醋酸钠 0 / 2.5 / 5 / 10 mM，各配等渗 NaCl 对照孔", "§5.6", TRACK_TOLERANCE),
    ],
    2: [
        Task("全量换液并重铺：400×g 5 min，重悬计数，调 3×10⁵/mL 回铺", "§5.3", TRACK_MAIN),
        Task("按组重新加药至工作浓度，不做残留折算", "§4.2 §5.3", TRACK_MAIN),
    ],
    3: [
        Task("耐受板计数并测活率", "§5.6", TRACK_TOLERANCE),
        Task("取活率降 <10% 且扩增降 <20% 的最高浓度（上限 5 mM），记下供第二轮用", "§5.6", TRACK_TOLERANCE),
    ],
    4: [
        Task("取样 → 立即流式染色", "§5.5 §6", TRACK_MAIN),
        Task("洗脱：PBS 洗 2 次，400×g 5 min", "§5.4", TRACK_MAIN),
        Task("对照组与撤药组重悬于仅含 IL-2 培养基；持续组重悬后加回 TGF-β 10 ng/mL", "§5.4", TRACK_MAIN),
        Task("计数，调回 3×10⁵/mL，1 mL/孔回铺", "§5.4", TRACK_MAIN),
    ],
    6: [
        Task("全量换液并重铺（同 Day 2）", "§5.3", TRACK_MAIN),
        Task("按组重新加药", "§4.2 §5.3", TRACK_MAIN),
    ],
    8: [
        Task("终点取样 → 立即流式染色", "§5.2 §5.5 §6", TRACK_MAIN),
        Task("核心六项：活死染料 / CD56 / GZMB / Perforin / CD103 / Ki67", "§6.2", TRACK_MAIN),
        Task("核对判定：撤药组 CD103 应回落至对照水平", "§8.1", TRACK_MAIN),
        Task("核对判定：撤药组 GZMB / Perforin 较对照低 ≥25%", "§8.1", TRACK_MAIN),
    ],
}

# 每次取样日都要做的记录动作（§9）
_SAMPLING_DAYS = (0, 4, 8)
_RECORD_TASKS = [
    Task("当天填表：相对天/组别/孔号/活细胞/死细胞/活率/传代数/操作者", "§9", TRACK_RECORD),
    Task("拍摄全板照片一张", "§9", TRACK_RECORD),
    Task(".fcs 按「日期_批次_组别」存入 raw/，只增不改", "§9", TRACK_RECORD),
]

# 开始前的准备。SOP 未指定日期，按 §2 §4.1 §6.1 汇总，自行安排。
PREP_TASKS = [
    Task("支原体检测", "§2.1"),
    Task("确认 FCS 单一批号，全程不换", "§2.1"),
    Task("复苏细胞，确保 Day 0 能达到 3×10⁵/mL；本轮独立重复 2 次需不同批次", "§5.1"),
    Task("TGF-β1 储液以 4 mM HCl + 0.1% BSA 配制并小体积分装 −80 °C", "§4.1"),
    Task("抗体 5 点两倍梯度滴定；关键标志阳性阈值以 FMO 确定", "§6.1"),
    Task("定下并记录板型（Day 4 用量决定是否改 12 / 6 孔板）", "§7.8"),
]

# 第二轮不排日期：它排在第一轮之后，且浓度依赖 §5.6 耐受板与 §7.2 浓度板。
ROUND2_GATES = [
    "先行 §7.2 浓度板（3 天）：醋酸钠取「能使全局 H3K27ac 升高的最低浓度」",
    "BMS-303141 取「H3K27ac 下降且活率降 <15%」的浓度；TSA 100 nM 4 h 作阳性对照",
    "第 2 组较第 1 组 GZMB / Perforin 须下降，否则返回第一轮，本轮无从判读",
    "承重比较是第 4 组（前段 Day 0–2）对第 5 组（后段 Day 2–4），暴露时长相同只差位置",
]


def next_weekday(after: _dt.date, weekday: int) -> _dt.date:
    """返回 ``after`` 之后（不含当天）的第一个指定星期几。0=周一。"""
    delta = (weekday - after.weekday() - 1) % 7 + 1
    return after + _dt.timedelta(days=delta)


def build_round1(day0: _dt.date) -> list[DayPlan]:
    """按 Day 0 生成第一轮全部有操作的天。"""
    plans: list[DayPlan] = []
    for rel_day in sorted(_ROUND1):
        tasks = list(_ROUND1[rel_day])
        if rel_day in _SAMPLING_DAYS:
            tasks.extend(_RECORD_TASKS)
        plans.append(DayPlan(rel_day=rel_day, date=day0 + _dt.timedelta(days=rel_day), tasks=tasks))
    return plans


def find_today(plans: list[DayPlan], today: _dt.date) -> DayPlan | None:
    for plan in plans:
        if plan.date == today:
            return plan
    return None


def find_next(plans: list[DayPlan], today: _dt.date) -> DayPlan | None:
    for plan in plans:
        if plan.date > today:
            return plan
    return None


def _date_label(day: _dt.date) -> str:
    return f"{day:%m-%d} 周{WEEKDAY_CN[day.weekday()]}"


# 空档达到这么多天才算「空得久」，值得提前预告
LONG_GAP = 2
# 提前几天开始预告
HEADS_UP_DAYS = 2


@dataclass
class Flow:
    """便签顶部那块随时间流动的面板所需的全部信息。"""

    today: _dt.date
    today_plan: DayPlan | None
    tomorrow_plan: DayPlan | None
    next_plan: DayPlan | None  # 严格晚于今天的下一个操作日
    days_until_next: int | None
    gap_len: int  # 今天所在这段空档有多少个连续无操作天

    @property
    def tomorrow(self) -> _dt.date:
        return self.today + _dt.timedelta(days=1)

    @property
    def heads_up(self) -> bool:
        """空档长、且离下一个操作日只剩 1–2 天时提前打招呼。

        空档只有一天时不提示 —— 「明天」那行已经把话说完了。
        """
        if self.today_plan is not None or self.days_until_next is None:
            return False
        return self.gap_len >= LONG_GAP and self.days_until_next <= HEADS_UP_DAYS


def _gap_length(plans: list[DayPlan], today: _dt.date, cap: int = 30) -> int:
    """今天所处空档的长度：往回数连续多少个无操作日，数到 ``cap`` 为止。

    只拿来和 ``LONG_GAP`` 比大小，所以封顶不影响判断；Day 0 之前空档无限长，
    封顶也省得一直往回走。
    """
    if find_today(plans, today) is not None:
        return 0
    operating = {p.date for p in plans}
    length, cursor = 1, today - _dt.timedelta(days=1)
    while length < cap and cursor not in operating:
        length += 1
        cursor -= _dt.timedelta(days=1)
    return length


def flow(plans: list[DayPlan], today: _dt.date) -> Flow:
    nxt = find_next(plans, today)
    return Flow(
        today=today,
        today_plan=find_today(plans, today),
        tomorrow_plan=find_today(plans, today + _dt.timedelta(days=1)),
        next_plan=nxt,
        days_until_next=(nxt.date - today).days if nxt else None,
        gap_len=_gap_length(plans, today),
    )


def today_line(f: Flow) -> str:
    head = f"今天 {_date_label(f.today)}"
    if f.today_plan is None:
        return f"{head} · 无操作"
    return f"{head} · {f.today_plan.label} · {len(f.today_plan.tasks)} 项"


def tomorrow_line(f: Flow) -> str:
    head = f"明天 {_date_label(f.tomorrow)}"
    if f.tomorrow_plan is not None:
        return f"{head} · {f.tomorrow_plan.label} · {len(f.tomorrow_plan.tasks)} 项"
    if f.next_plan is None:
        return f"{head} · 无操作 · 第一轮已结束"
    # 明天是等着的，就把要等到几号说清楚
    return (
        f"{head} · 等待 · 下一次 {f.next_plan.label} 在 "
        f"{_date_label(f.next_plan.date)}，还有 {f.days_until_next} 天"
    )


@dataclass(frozen=True)
class Notification:
    title: str
    subtitle: str
    body: str


def notification_for(f: Flow) -> Notification | None:
    """当天该不该弹通知，弹什么。没什么可说的就返回 None，不制造噪音。"""
    if f.today_plan is not None:
        first = f.today_plan.tasks[0].text if f.today_plan.tasks else ""
        return Notification(
            title=f"TGF-β SOP · {f.today_plan.label}",
            subtitle=f"今天 {_date_label(f.today)} · {len(f.today_plan.tasks)} 项",
            body=first,
        )
    if f.heads_up and f.next_plan is not None:
        when = "明天" if f.days_until_next == 1 else f"还有 {f.days_until_next} 天"
        first = f.next_plan.tasks[0].text if f.next_plan.tasks else ""
        return Notification(
            title="TGF-β SOP · 预告",
            subtitle=f"{when}就是 {f.next_plan.label}（{_date_label(f.next_plan.date)}）",
            body=first,
        )
    return None


def heads_up_line(f: Flow) -> str | None:
    """空档长时的提前预告；不该提示就返回 None。"""
    if not f.heads_up or f.next_plan is None:
        return None
    when = "明天" if f.days_until_next == 1 else f"还有 {f.days_until_next} 天"
    first = f.next_plan.tasks[0].text if f.next_plan.tasks else ""
    return f"⏰ {when}就是 {f.next_plan.label}（{_date_label(f.next_plan.date)}）：{first}"


def render_text(day0: _dt.date, today: _dt.date | None = None) -> str:
    """纯文本版日程，供打印或贴进实验记录本。"""
    today = today or _dt.date.today()
    plans = build_round1(day0)
    out = [
        SOP_TITLE,
        f"第一轮：撤药模型建立　Day 0 = {day0:%Y-%m-%d}（周{WEEKDAY_CN[day0.weekday()]}）",
        "",
        "开始前（SOP 未指定日期，自行安排）",
    ]
    out += [f"  · {t.text}　{t.ref}" for t in PREP_TASKS]
    out.append("")
    for plan in plans:
        out.append(f"{plan.label}　{plan.date:%Y-%m-%d} 周{WEEKDAY_CN[plan.date.weekday()]}")
        for task in plan.tasks:
            tag = "" if task.track == TRACK_MAIN else f"[{TRACK_LABEL[task.track]}] "
            out.append(f"  · {tag}{task.text}　{task.ref}")
        out.append("")
    out.append("第二轮：窄通道检验　——　不排日期，条件性启动")
    out += [f"  · {g}" for g in ROUND2_GATES]
    return "\n".join(out)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="打印 SOP 第一轮日程")
    parser.add_argument("--day0", help="Day 0 日期 YYYY-MM-DD；缺省为下一个周二")
    args = parser.parse_args()
    if args.day0:
        start = _dt.date.fromisoformat(args.day0)
    else:
        start = next_weekday(_dt.date.today(), 1)
    print(render_text(start))
