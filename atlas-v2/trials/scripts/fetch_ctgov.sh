#!/bin/bash
# usage: fetch_ctgov.sh <outprefix> <cond query>
out=$1; q=$2; tok=""; i=0
while :; do
  url="https://clinicaltrials.gov/api/v2/studies?query.cond=$(python3 -c 'import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))' "$q")&filter.advanced=AREA%5BStudyType%5DINTERVENTIONAL&pageSize=200&countTotal=true"
  [ -n "$tok" ] && url="$url&pageToken=$tok"
  curl -sS "$url" -o "${out}_$i.json" || exit 1
  tok=$(python3 -I -c 'import json,sys;print(json.load(open(sys.argv[1])).get("nextPageToken",""))' "${out}_$i.json")
  i=$((i+1)); [ -z "$tok" ] && break
done
echo "$out pages=$i"
