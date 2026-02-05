from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from mutme.utils import (
    CommandExecutionError,
    run_command,
)


def py_cmd(code: str) -> list[str]:
    """Portable helper: run a tiny Python snippet."""
    return [sys.executable, "-c", code]


def test_success_captures_stdout_and_stderr():
    result = run_command(py_cmd('import sys; print("out"); print("err", file=sys.stderr)'))
    assert result.returncode == 0
    assert result.stdout.strip() == "out"
    assert result.stderr.strip() == "err"
    assert result.cwd is None
    assert isinstance(result.args, tuple)
    assert result.args[0] == sys.executable


def test_capture_output_false_returns_empty_strings():
    result = run_command(py_cmd('print("hello")'), capture_output=False)
    assert result.returncode == 0
    assert result.stdout == ""
    assert result.stderr == ""


def test_check_true_raises_rich_error_on_nonzero_exit():
    with pytest.raises(CommandExecutionError) as excinfo:
        run_command(py_cmd("raise SystemExit(2)"), check=True)

    err = excinfo.value
    assert "exit code 2" in str(err)
    assert err.result is not None
    assert err.result.returncode == 2
    assert err.result.args[0] == sys.executable


def test_check_false_returns_result_on_nonzero_exit():
    result = run_command(py_cmd("raise SystemExit(7)"), check=False)
    assert result.returncode == 7


def test_command_str_splits_when_shell_false():
    # This uses your shlex.split path, but still runs portably.
    cmd = f'{sys.executable} -c "print(123)"'
    result = run_command(cmd, shell=False)
    assert result.returncode == 0
    assert result.stdout.strip() == "123"
    # args should be the split tuple; first element is the interpreter path
    assert result.args[0] == sys.executable


def test_command_sequence_requires_all_str():
    with pytest.raises(TypeError, match="sequence of str"):
        run_command([sys.executable, 123])  # type: ignore[list-item]


def test_command_invalid_type_rejected():
    with pytest.raises(TypeError, match="str or a sequence of str"):
        run_command(123)  # type: ignore[arg-type]


def test_cwd_validation_nonexistent(tmp_path: Path):
    missing = tmp_path / "nope"
    with pytest.raises(ValueError, match="cwd does not exist"):
        run_command(py_cmd("print('x')"), cwd=missing)


def test_cwd_validation_not_a_directory(tmp_path: Path):
    file_path = tmp_path / "file.txt"
    file_path.write_text("hi", encoding="utf-8")

    with pytest.raises(ValueError, match="cwd is not a directory"):
        run_command(py_cmd("print('x')"), cwd=file_path)


def test_cwd_is_applied(tmp_path: Path):
    # Print current working directory from inside the child process
    result = run_command(py_cmd("import os; print(os.getcwd())"), cwd=tmp_path)
    assert result.returncode == 0
    assert Path(result.stdout.strip()) == tmp_path
    assert result.cwd == str(tmp_path)


def test_env_overlay_is_applied():
    result = run_command(
        py_cmd('import os; print(os.getenv("MY_TEST_VAR", ""))'),
        env={"MY_TEST_VAR": "hello"},
    )
    assert result.returncode == 0
    assert result.stdout.strip() == "hello"


def test_stdin_text_mode_roundtrip():
    result = run_command(
        py_cmd("import sys; data=sys.stdin.read(); print(data, end='')"),
        stdin="hello stdin",
        text=True,
    )
    assert result.returncode == 0
    assert result.stdout == "hello stdin"


def test_stdin_bytes_with_text_true_rejected():
    with pytest.raises(TypeError, match="stdin is bytes but text=True"):
        run_command(py_cmd("print('x')"), stdin=b"nope", text=True)


def test_stdin_str_with_text_false_rejected():
    with pytest.raises(TypeError, match="stdin is str but text=False"):
        run_command(py_cmd("print('x')"), stdin="nope", text=False)


def test_file_not_found_is_wrapped(monkeypatch: pytest.MonkeyPatch):
    def boom(*args, **kwargs):
        raise FileNotFoundError("no such file")

    monkeypatch.setattr(subprocess, "run", boom)

    with pytest.raises(CommandExecutionError, match="Command not found"):
        run_command(["definitely-not-a-real-command-xyz"])


def test_timeout_is_wrapped_and_includes_partial_result(
    monkeypatch: pytest.MonkeyPatch,
):
    # Simulate subprocess.TimeoutExpired with partial stdout/stderr
    def boom(*args, **kwargs):
        raise subprocess.TimeoutExpired(
            cmd=kwargs.get("args", args[0]),
            timeout=0.01,
            output=b"partial out",
            stderr=b"partial err",
        )

    monkeypatch.setattr(subprocess, "run", boom)

    with pytest.raises(CommandExecutionError) as excinfo:
        run_command(["something"], timeout=0.01)

    err = excinfo.value
    assert "timed out" in str(err)
    assert err.result is not None
    assert err.result.returncode == -1
    assert "partial out" in err.result.stdout
    assert "partial err" in err.result.stderr


def test_generic_oserror_is_wrapped(monkeypatch: pytest.MonkeyPatch):
    def boom(*args, **kwargs):
        raise OSError("weird exec failure")

    monkeypatch.setattr(subprocess, "run", boom)

    with pytest.raises(CommandExecutionError, match="Failed to execute command"):
        run_command(["something"])
