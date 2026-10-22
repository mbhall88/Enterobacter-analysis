"""Run one command, keeping argv, timing, exit status and output in logs/NAME.*."""
import datetime as dt
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

LOGS = Path(__file__).resolve().parents[1] / "logs"


def main():
    name, argv = sys.argv[1], sys.argv[2:]
    record = LOGS / f"{name}.json"
    if record.exists():
        raise SystemExit(f"refusing to overwrite {record}")
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    clock = time.monotonic()
    with (LOGS / f"{name}.stdout").open("x") as out, (LOGS / f"{name}.stderr").open("x") as err:
        code = subprocess.run(argv, stdout=out, stderr=err, check=False).returncode
    with record.open("x") as handle:
        json.dump(dict(argv=argv, cwd=str(Path.cwd()), started_utc=started, elapsed_seconds=time.monotonic() - clock,
                       exit_code=code, child_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),
                  handle, indent=2)
        handle.write("\n")
    print((LOGS / f"{name}.stdout").read_text()[-4000:])
    print((LOGS / f"{name}.stderr").read_text()[-4000:], file=sys.stderr)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
