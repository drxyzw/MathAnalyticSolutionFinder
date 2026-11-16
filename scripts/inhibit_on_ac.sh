#!/bin/bash
# Hold a sleep/idle inhibitor only while the laptop is on external power, for as long as <pid> is alive.
# On battery the inhibitor is released, so the normal sleep settings apply (the run is suspended, not killed).
# Usage: setsid nohup scripts/inhibit_on_ac.sh <pid> > <log> 2>&1 < /dev/null &
PID=$1
[ -z "$PID" ] && { echo "usage: $0 <pid>"; exit 1; }
INH=

on_ac() {
	# "Discharging" means running on battery; Charging / Full / Not charging mean external power.
	local s
	s=$(cat /sys/class/power_supply/BAT*/status 2>/dev/null | head -1)
	[ -z "$s" ] && return 0   # no battery: treat as AC
	[ "$s" != "Discharging" ]
}

release() {
	if [ -n "$INH" ]; then
		pkill -P "$INH" 2>/dev/null
		kill "$INH" 2>/dev/null
		INH=
	fi
}
trap release EXIT

while kill -0 "$PID" 2>/dev/null; do
	if on_ac; then
		if [ -z "$INH" ] || ! kill -0 "$INH" 2>/dev/null; then
			systemd-inhibit --what=sleep:idle --who="zeta PSLQ" --why="production run (AC power only)" --mode=block sleep infinity &
			INH=$!
			echo "$(date) AC power: inhibitor ON"
		fi
	else
		if [ -n "$INH" ]; then
			release
			echo "$(date) battery: inhibitor OFF"
		fi
	fi
	sleep 30
done
echo "$(date) pid $PID finished; exiting"
