"""把提醒送进 macOS 通知中心。

用系统自带的 osascript ``display notification``，不装第三方包。
其他平台上 :func:`send` 直接返回 False，便签本身照常用。
"""

from __future__ import annotations

import shutil
import subprocess
import sys

TIMEOUT = 10  # osascript 卡住时别把便签一起拖死


def available() -> bool:
    """当前系统能不能发通知。"""
    return sys.platform == "darwin" and shutil.which("osascript") is not None


def applescript_string(text: str) -> str:
    """把 Python 字符串包成 AppleScript 字面量。

    AppleScript 的双引号字符串里，反斜杠和双引号都要转义；换行会截断字面量，
    折成空格。反斜杠必须先换，否则会把后面补的反斜杠再转义一遍。
    """
    escaped = text.replace("\\", "\\\\").replace('"', '\\"')
    escaped = escaped.replace("\r\n", " ").replace("\n", " ").replace("\r", " ")
    return f'"{escaped}"'


def build_script(title: str, subtitle: str, body: str) -> str:
    return (
        f"display notification {applescript_string(body)} "
        f"with title {applescript_string(title)} "
        f"subtitle {applescript_string(subtitle)}"
    )


def send(title: str, subtitle: str, body: str) -> bool:
    """发一条通知。成功返回 True；发不出去返回 False，绝不抛异常。

    通知被系统关掉时 osascript 往往照样返回 0，所以 True 只代表命令跑通了，
    不代表用户真看见了。
    """
    if not available():
        return False
    try:
        done = subprocess.run(
            ["osascript", "-e", build_script(title, subtitle, body)],
            capture_output=True,
            timeout=TIMEOUT,
            check=False,
        )
        return done.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False
