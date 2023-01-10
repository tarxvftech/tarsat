
#!/bin/bash
epoch="$(cat $1 |head -n 2 |tail -n 1|awk '{print $4;}'|sed 's/\..*//')"
echo $epoch
if [[ ! -z "$epoch" ]]; then
	cp $1 historical/$epoch.txt
fi
