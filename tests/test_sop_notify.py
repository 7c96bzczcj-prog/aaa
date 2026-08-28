"""Tests for the macOS notification path (sop_memo/notify.py + memo.maybe_notify).

投递本身要 macOS，这里跑在别的系统上，所以锁的是能锁的部分：AppleScript 字面量
的转义、脚本拼装、发不出去时不抛异常，以及「每天最多一条、到点才发、没内容不发」
这套判断。
"""

from __future__ import annotations

import datetime as dt
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "sop_memo"))

import notify  # noqa: E402
from schedule import build_round1, flow, notification_for  # noqa: E402

DAY0 = dt.date(2026, 9, 1)


# ------------------------------------------------- AppleScript 字面量转义
def test_plain_text_is_just_quoted():
    assert notify.applescript_string("铺板") == '"铺板"'


def test_double_quotes_are_escaped():
    assert notify.applescript_string('说"你好"') == '"说\\"你好\\""'


def test_backslashes_are_escaped_before_quotes():
    """反斜杠必须先转义，否则会把补出来的反斜杠再转义一遍。"""
    assert notify.applescript_string(r"a\b") == r'"a\\b"'
    assert notify.applescript_string('\\"') == '"\\\\\\""'


def test_newlines_are_folded_to_spaces():
    """换行会截断 AppleScript 字面量，得折平。"""
    assert notify.applescript_string("一\n二\r\n三\r四") == '"一 二 三 四"'


def test_build_script_puts_all_three_fields_in():
    script = notify.build_script("标题", "副标题", "正文")
    assert script == (
        'display notification "正文" with title "标题" subtitle "副标题"'
    )


def test_build_script_survives_a_quote_in_the_body():
    """带引号的正文不能把脚本拼断。"""
    script = notify.build_script("t", "s", '他说"做完了"')
    assert script.count('"') % 2 == 0
    assert '\\"做完了\\"' in script


# ------------------------------------------------------------ send 的容错
def test_send_is_a_noop_off_macos(monkeypatch):
    monkeypatch.setattr(sys, "platform", "linux")
    assert notify.available() is False
    assert notify.send("t", "s", "b") is False


def test_send_reports_failure_instead_of_raising(monkeypatch):
    """osascript 缺失或超时都不能把便签带崩。"""
    monkeypatch.setattr(notify, "available", lambda: True)

    def boom(*_a, **_kw):
        raise OSError("no osascript")

    monkeypatch.setattr(notify.subprocess, "run", boom)
    assert notify.send("t", "s", "b") is False

    def slow(*_a, **_kw):
        raise subprocess.TimeoutExpired(cmd="osascript", timeout=10)

    monkeypatch.setattr(notify.subprocess, "run", slow)
    assert notify.send("t", "s", "b") is False


def test_send_passes_argv_without_a_shell(monkeypatch):
    """走 argv，不经 shell —— 文案里的字符不该被当成 shell 语法。"""
    seen = {}
    monkeypatch.setattr(notify, "available", lambda: True)

    def fake(cmd, **kwargs):
        seen["cmd"], seen["kwargs"] = cmd, kwargs
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(notify.subprocess, "run", fake)
    assert notify.send("标题", "副标题", "正文") is True
    assert seen["cmd"][:2] == ["osascript", "-e"]
    assert len(seen["cmd"]) == 3
    assert "shell" not in seen["kwargs"]


def test_send_reports_a_nonzero_exit_as_failure(monkeypatch):
    monkeypatch.setattr(notify, "available", lambda: True)
    monkeypatch.setattr(
        notify.subprocess, "run",
        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1),
    )
    assert notify.send("t", "s", "b") is False


# ------------------------------------------------------- 通知内容怎么定
def _note(today: dt.date):
    return notification_for(flow(build_round1(DAY0), today))


def test_an_operating_day_notifies_with_the_day_and_count():
    note = _note(DAY0)
    assert note.title == "TGF-β SOP · Day 0"
    assert note.subtitle == "今天 09-01 周二 · 8 项"
    assert "铺 24 孔板" in note.body


def test_the_heads_up_day_notifies_with_the_countdown():
    note = _note(dt.date(2026, 8, 30))  # 离 Day 0 两天
    assert note.title == "TGF-β SOP · 预告"
    assert note.subtitle == "还有 2 天就是 Day 0（09-01 周二）"


def test_a_quiet_day_notifies_nothing():
    """离得远、或空档只有一天，都不打扰。"""
    assert _note(dt.date(2026, 8, 28)) is None  # 离 Day 0 还有 4 天
    assert _note(DAY0 + dt.timedelta(days=1)) is None  # Day 1，明天那行已说明
    assert _note(DAY0 + dt.timedelta(days=20)) is None  # 第一轮结束之后


# --------------------------------------------- 每天一条 / 到点才发 的判断
class FakeMemo:
    """只取 Memo.maybe_notify 需要的那几样，避开 tkinter。"""

    def __init__(self, notify_at=dt.time(9, 0)):
        self.notify_at = notify_at
        self.plans = build_round1(DAY0)
        self.state = {"checked": {}}
        self.sent = []

    maybe_notify = None  # 由 fixture 绑上真方法


SENT: list[tuple[str, str, str]] = []


@pytest.fixture
def memo_cls(monkeypatch):
    import memo as memo_module

    SENT.clear()

    def record(title, subtitle, body):
        SENT.append((title, subtitle, body))
        return True

    monkeypatch.setattr(memo_module, "save_state", lambda _s: None)
    monkeypatch.setattr(memo_module.notify, "send", record)
    FakeMemo.maybe_notify = memo_module.Memo.maybe_notify
    return FakeMemo


def test_the_day_content_is_what_actually_reaches_the_notifier(memo_cls):
    memo_cls().maybe_notify(dt.datetime(2026, 9, 5, 9, 0))  # Day 4
    assert len(SENT) == 1
    title, subtitle, body = SENT[0]
    assert title == "TGF-β SOP · Day 4"
    assert subtitle == "今天 09-05 周六 · 7 项"
    assert body == "取样 → 立即流式染色"


def test_nothing_fires_before_the_appointed_time(memo_cls):
    m = memo_cls()
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 8, 59)) is False
    assert "last_notified" not in m.state


def test_it_fires_once_the_time_has_come(memo_cls):
    m = memo_cls()
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 9, 0)) is True
    assert m.state["last_notified"] == "2026-09-01"


def test_it_fires_only_once_a_day(memo_cls):
    m = memo_cls()
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 9, 0)) is True
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 9, 5)) is False
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 23, 59)) is False


def test_a_new_day_fires_again(memo_cls):
    m = memo_cls()
    m.maybe_notify(dt.datetime(2026, 9, 1, 9, 0))
    assert m.maybe_notify(dt.datetime(2026, 9, 3, 9, 0)) is True  # Day 2
    assert m.state["last_notified"] == "2026-09-03"


def test_a_late_start_still_gets_the_days_notice(memo_cls):
    """中午才开机，当天那条要补上，不是跳过。"""
    m = memo_cls()
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 14, 30)) is True


def test_a_quiet_day_marks_itself_done_without_sending(memo_cls):
    """没内容也记一笔，免得每 5 分钟重算。"""
    m = memo_cls()
    assert m.maybe_notify(dt.datetime(2026, 9, 2, 9, 0)) is False  # Day 1 空
    assert m.state["last_notified"] == "2026-09-02"


def test_no_notify_switches_it_off_entirely(memo_cls):
    m = memo_cls(notify_at=None)
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 9, 0)) is False
    assert "last_notified" not in m.state


def test_a_custom_time_is_respected(memo_cls):
    m = memo_cls(notify_at=dt.time(7, 30))
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 7, 29)) is False
    assert m.maybe_notify(dt.datetime(2026, 9, 1, 7, 30)) is True


# ------------------------------------------------------------ --notify-at
def test_notify_at_parses_and_rejects_junk():
    import memo as memo_module

    assert memo_module.parse_args(["--notify-at", "07:30"]).notify_at == dt.time(7, 30)
    with pytest.raises(SystemExit):
        memo_module.parse_args(["--notify-at", "半夜"])
