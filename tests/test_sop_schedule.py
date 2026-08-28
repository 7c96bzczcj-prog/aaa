"""Tests for the SOP memo schedule (sop_memo/schedule.py).

相对天来自 SOP §5，不可改；这里锁住的是相对天 → 日历日期的换算，
以及取样日必须带上 §9 的记录动作。
"""

from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "sop_memo"))

from schedule import (  # noqa: E402
    TRACK_RECORD,
    TRACK_TOLERANCE,
    build_round1,
    find_next,
    find_today,
    flow,
    heads_up_line,
    next_weekday,
    render_text,
    today_line,
    tomorrow_line,
)

FRIDAY = dt.date(2026, 8, 28)  # 交代任务当天
DAY0 = dt.date(2026, 9, 1)  # 下周二


# ------------------------------------------------------------- next_weekday
def test_next_tuesday_from_the_friday_the_task_was_given():
    assert next_weekday(FRIDAY, 1) == DAY0
    assert DAY0.strftime("%A") == "Tuesday"


def test_next_weekday_is_strictly_after_today():
    """今天就是周二时，要给出下一个周二，而不是今天。"""
    assert next_weekday(DAY0, 1) == DAY0 + dt.timedelta(days=7)


# -------------------------------------------------------------- build_round1
def test_operating_days_match_the_sop():
    """SOP §5.2 排 Day 0/2/4/6/8，§5.6 的耐受板在 Day 3 加一天。"""
    plans = build_round1(DAY0)
    assert [p.rel_day for p in plans] == [0, 2, 3, 4, 6, 8]


def test_relative_days_map_onto_the_right_dates():
    plans = {p.rel_day: p.date for p in build_round1(DAY0)}
    assert plans[0] == dt.date(2026, 9, 1)
    assert plans[2] == dt.date(2026, 9, 3)
    assert plans[3] == dt.date(2026, 9, 4)
    assert plans[4] == dt.date(2026, 9, 5)
    assert plans[6] == dt.date(2026, 9, 7)
    assert plans[8] == dt.date(2026, 9, 9)


def test_sampling_days_carry_the_record_actions():
    """§5.5 的取样日是 Day 0/4/8，每个都要带 §9 的填表、拍照、存档。"""
    plans = {p.rel_day: p for p in build_round1(DAY0)}
    for rel_day in (0, 4, 8):
        tracks = [t.track for t in plans[rel_day].tasks]
        assert tracks.count(TRACK_RECORD) == 3, f"Day {rel_day} 少了记录动作"
    for rel_day in (2, 6):
        assert TRACK_RECORD not in [t.track for t in plans[rel_day].tasks]


def test_tolerance_plate_opens_on_day0_and_reads_out_on_day3():
    plans = {p.rel_day: p for p in build_round1(DAY0)}
    assert any(t.track == TRACK_TOLERANCE for t in plans[0].tasks)
    assert all(t.track == TRACK_TOLERANCE for t in plans[3].tasks)


def test_day4_washout_keeps_all_four_steps():
    """§5.4 洗脱是四步：转管洗涤、按组重悬、计数回铺，外加取样。"""
    day4 = {p.rel_day: p for p in build_round1(DAY0)}[4]
    texts = " ".join(t.text for t in day4.tasks)
    assert "PBS 洗 2 次" in texts
    assert "仅含 IL-2" in texts
    assert "3×10⁵/mL" in texts


# ---------------------------------------------------- 顶部面板：今天 / 明天
def _panel(today: dt.date):
    f = flow(build_round1(DAY0), today)
    return today_line(f), tomorrow_line(f), heads_up_line(f)


def test_today_line_names_the_day_and_counts_the_work():
    today, _, _ = _panel(DAY0)
    assert today == "今天 09-01 周二 · Day 0 · 8 项"


def test_today_line_says_so_when_there_is_nothing_to_do():
    today, _, _ = _panel(DAY0 + dt.timedelta(days=1))  # Day 1 空
    assert today == "今天 09-02 周三 · 无操作"


def test_tomorrow_line_names_tomorrows_work():
    _, tomorrow, _ = _panel(DAY0 + dt.timedelta(days=1))  # 明天是 Day 2
    assert tomorrow == "明天 09-03 周四 · Day 2 · 2 项"


def test_tomorrow_line_spells_out_the_date_when_tomorrow_is_a_wait():
    """用户要求：明天如果是等待，得写明等到几号。"""
    _, tomorrow, _ = _panel(FRIDAY)  # 08-29 周六仍在等 Day 0
    assert tomorrow == "明天 08-29 周六 · 等待 · 下一次 Day 0 在 09-01 周二，还有 4 天"


def test_tomorrow_line_reports_the_round_as_finished():
    _, tomorrow, _ = _panel(DAY0 + dt.timedelta(days=9))
    assert tomorrow.endswith("第一轮已结束")


# ------------------------------------------------------ 空档长时的提前预告
def test_no_heads_up_while_the_gap_is_still_far_off():
    """08-28 离 Day 0 还有 4 天，不用现在就催。"""
    assert _panel(FRIDAY)[2] is None


def test_heads_up_fires_two_days_out_after_a_long_gap():
    *_, ahead = _panel(dt.date(2026, 8, 30))  # 离 Day 0 两天
    assert ahead is not None
    assert "还有 2 天就是 Day 0（09-01 周二）" in ahead
    assert "铺 24 孔板" in ahead  # 带上那天头一件事，不是光报个日子


def test_heads_up_fires_one_day_out_too():
    *_, ahead = _panel(dt.date(2026, 8, 31))
    assert ahead is not None and ahead.startswith("⏰ 明天就是 Day 0")


def test_no_heads_up_on_an_operating_day():
    assert _panel(DAY0)[2] is None


def test_a_single_idle_day_gets_no_heads_up():
    """Day 1 是唯一的空档日，「明天」那行已经说清楚了，不再重复催。"""
    plans = build_round1(DAY0)
    f = flow(plans, DAY0 + dt.timedelta(days=1))
    assert f.gap_len == 1
    assert f.days_until_next == 1
    assert heads_up_line(f) is None


def test_flow_reports_the_gap_and_the_distance():
    f = flow(build_round1(DAY0), FRIDAY)
    assert f.today_plan is None
    assert f.next_plan.rel_day == 0
    assert f.days_until_next == 4
    assert f.gap_len >= 2  # Day 0 之前是长空档


# ------------------------------------------------------------ find_today/next
def test_find_today_and_find_next():
    plans = build_round1(DAY0)
    assert find_today(plans, dt.date(2026, 9, 3)).rel_day == 2
    assert find_today(plans, dt.date(2026, 9, 2)) is None
    assert find_next(plans, dt.date(2026, 9, 2)).rel_day == 2
    assert find_next(plans, dt.date(2026, 9, 9)) is None


# ------------------------------------------------------------------ rendering
def test_text_render_carries_dates_prep_and_the_round2_gate():
    out = render_text(DAY0, today=FRIDAY)
    assert "Day 0 = 2026-09-01（周二）" in out
    assert "Day 4　2026-09-05 周六" in out
    assert "支原体检测" in out  # 开始前
    assert "第二轮" in out and "条件性启动" in out
