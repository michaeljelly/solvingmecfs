#!/bin/bash
# usage: fetch_epmc.sh <outprefix> <query>
out=$1; q=$2; cm="*"; i=0
while :; do
  curl -sS -G "https://www.ebi.ac.uk/europepmc/webservices/rest/search" --data-urlencode "query=$q" --data-urlencode "format=json" --data-urlencode "pageSize=1000" --data-urlencode "resultType=core" --data-urlencode "cursorMark=$cm" -o "${out}_$i.json" || exit 1
  n=$(python3 -I -c 'import json,sys;d=json.load(open(sys.argv[1]));print(len(d["resultList"]["result"]))' "${out}_$i.json")
  ncm=$(python3 -I -c 'import json,sys;d=json.load(open(sys.argv[1]));print(d.get("nextCursorMark",""))' "${out}_$i.json")
  i=$((i+1))
  [ "$n" -lt 1000 ] && break; [ "$ncm" = "$cm" ] && break; cm=$ncm
done
echo "$out pages=$i"
