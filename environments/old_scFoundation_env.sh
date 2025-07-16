mamba create --prefix /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env python=3.10 -y
conda activate /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env
pip install numpy
pip install pandas
pip install scipy
pip install torch==2.5.0+cu121 --index-url https://download.pytorch.org/whl/cu121
pip install einops
pip install scanpy
pip install local_attention
#Download model from this URL: https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FPublicSharedfiles%2FShared%20Documents%2FPublic%20Shared%20files&p=true&ga=1
mkdir -p /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models
# Use the following curl request 
curl 'https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/_layouts/15/download.aspx?UniqueId=f19a9c0b%2Dffac%2D458c%2Dbe21%2D4c2db2f3e17b' \
  -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:139.0) Gecko/20100101 Firefox/139.0' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Referer: https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FPublicSharedfiles%2FShared%20Documents%2FPublic%20Shared%20files^&p=true^&ga=1' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'Sec-Fetch-Dest: iframe' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Connection: keep-alive' \
  -H 'Cookie: rtFa=YfXpo1sKB0qOwOSTCtTNGF6pSafxRvDhmLuh+pOTg8wmM2I5ODY1YTAtZjA2Ni00ODQ2LWJmNGQtMTEzZjRiZWYwNzU1IzEzMzkzMzQwNzA4OTIyMjUxMSM3MjJjYTRhMS1jMDdlLTkwMDAtMGI1Yi1lOGI0YjBlZDFiNzYjbmthbGF2cm8lNDBiaWRtYy5oYXJ2YXJkLmVkdSMxOTMzMTIjcDI4SXlmdEpQaGgydnlRbUFjdmZSR0xSRjU0I3l2ZWFaSzZ4TzVCbXJDUE9KQlE0OWoyRDE2axCpK03g714/joCUroXJ+ly9w6mR0/IqKTp7gRnb91zWK6Qg9P8ETeN5ITl8oIDW/azcHkO2y/lSVRwe9dSFVUW6PqGv8EVRS/fnjtsdX0mzAgjFLTqzO7O7f3aYYmtNWayPudr6YxB+0S06KZnHgj4q8whKV31vi9+RLjcTFlWyfi/6SFFHQx2bZkKKxkETreEEzYfi4WiJiXojy/Tzhxkm617UsG/yzNjXn6EJ+a1lBF5vfyB0I6vQ8rKqSWcj/bl1Pm/bRyAAF86BVGuEgEU6t3iQKzSXJMFDV8Gi7eTJE1q2IwSWifz/7Z+lwzvSYmF0oapFKtBX6LLcNVZXVcrZAAAA; SIMI=eyJzdCI6MTc0NjA0NDczMjgzNH0=; spo_abt=Mjg4MDAsLFsiYXBwX3JlcyJdLDAwM2Y4MWU5LTRjYTMtNjhmOC1iYTQ3LTg3Yzg2NzcxNTQyNQ==; FedAuth=77u/PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0idXRmLTgiPz48U1A+VjE0LDBoLmZ8bWVtYmVyc2hpcHx1cm4lM2FzcG8lM2Fhbm9uIzUyOWY3NjQ1ZWFhOTY1ZjkxNzcyYmJhNWQ3M2E1ZjczMDdkNWY3ZTQzZTgyOTBkODA5NmM2ZmU2NTM5MTY0ZmYsMCMuZnxtZW1iZXJzaGlwfHVybiUzYXNwbyUzYWFub24jNTI5Zjc2NDVlYWE5NjVmOTE3NzJiYmE1ZDczYTVmNzMwN2Q1ZjdlNDNlODI5MGQ4MDk2YzZmZTY1MzkxNjRmZiwxMzM5NDY2MTMwNDAwMDAwMDAsMCwxMzM5NDc0NzQwNTA1ODAyNTksMC4wLjAuMCwyNTgsYjk1ZjUyMjgtNjM5Yy00ZmViLWEwMDQtNzcwYWRjMTQwNjdlLCwsNWM4MjYwMWEtZjI1ZS00YWQxLTg5NDktMmEwMjcyZTIwOWY3LDVjODI2MDFhLWYyNWUtNGFkMS04OTQ5LTJhMDI3MmUyMDlmNyxMWmpoMjVWL1gwV0ZsRUhVSVFuUEpnLDAsMCwwLCwsLDI2NTA0Njc3NDM5OTk5OTk5OTksMCwsLCwsLCwwLCwxOTU4OTcsbWFmVnZaa0RRMTl1WkdIVnhRYWpXNkZLMDBJLCxHb24wS2xHbEEzaVFidjBkcnU4RkJsTU9FdGJWV2tnYW02eW5rZjVQOTkya01pYmYzQk10cTFZcWJNbnN3SExrTVM1bHpIUzhySjM0enlXUG9ydVNYMjJ1TmpjNVoyLytQekZPcWVuUVVnT01YVG1TSUpUVmgvR2dIU2VnKzNlM1lUTjNPWno0Wk1xdDdoM0ZONEJiWTAyMEZQQ3JPSnMvQzlyOFVzMFhZMXVRODlRZUlSTWVjTkI0QWhMMFJlLzNHTS90ZnZDeDV1WnVoWWtiN3lET2g0c2JMRk5uS3FtZC9jWG5jbzFuZi9QcCtuV1FFNUtvRTdTSTc2Sm1VcVcyU1lPWW13TjVpRmlpL093b1hjSklwMWg3bFpUamdVcjFWdmJabzBST0hKMGw1Zmo2bWtXYzdvZWl4RVlNZUJxN09ZUTNQYitVR1pUaUQrY21JbWlKbWc9PTwvU1A+; FeatureOverrides_experiments=[]; MicrosoftApplicationsTelemetryDeviceId=a5ab09b6-1a7f-467e-8b93-10ae64080871; ai_session=QHl+appkL0i1W8CeLuMONI^|1750187406931^|1750188963111; MSFPC=GUID=fc9103de21f54b1fb816e508597e4fa5^&HASH=fc91^&LV=202502^&V=4^&LU=1740503289506; SPA_RT=' \
  -o /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models/models.ckpt

curl 'https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/_layouts/15/download.aspx?UniqueId=f7a65086%2Da179%2D405d%2Dbd99%2D397675621e83' \
  -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:139.0) Gecko/20100101 Firefox/139.0' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Referer: https://hopebio2020.sharepoint.com/sites/PublicSharedfiles/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FPublicSharedfiles%2FShared%20Documents%2FPublic%20Shared%20files^&p=true^&ga=1' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'Sec-Fetch-Dest: iframe' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Connection: keep-alive' \
  -H 'Cookie: rtFa=YfXpo1sKB0qOwOSTCtTNGF6pSafxRvDhmLuh+pOTg8wmM2I5ODY1YTAtZjA2Ni00ODQ2LWJmNGQtMTEzZjRiZWYwNzU1IzEzMzkzMzQwNzA4OTIyMjUxMSM3MjJjYTRhMS1jMDdlLTkwMDAtMGI1Yi1lOGI0YjBlZDFiNzYjbmthbGF2cm8lNDBiaWRtYy5oYXJ2YXJkLmVkdSMxOTMzMTIjcDI4SXlmdEpQaGgydnlRbUFjdmZSR0xSRjU0I3l2ZWFaSzZ4TzVCbXJDUE9KQlE0OWoyRDE2axCpK03g714/joCUroXJ+ly9w6mR0/IqKTp7gRnb91zWK6Qg9P8ETeN5ITl8oIDW/azcHkO2y/lSVRwe9dSFVUW6PqGv8EVRS/fnjtsdX0mzAgjFLTqzO7O7f3aYYmtNWayPudr6YxB+0S06KZnHgj4q8whKV31vi9+RLjcTFlWyfi/6SFFHQx2bZkKKxkETreEEzYfi4WiJiXojy/Tzhxkm617UsG/yzNjXn6EJ+a1lBF5vfyB0I6vQ8rKqSWcj/bl1Pm/bRyAAF86BVGuEgEU6t3iQKzSXJMFDV8Gi7eTJE1q2IwSWifz/7Z+lwzvSYmF0oapFKtBX6LLcNVZXVcrZAAAA; SIMI=eyJzdCI6MTc0NjA0NDczMjgzNH0=; spo_abt=Mjg4MDAsLFsiYXBwX3JlcyJdLDAwM2Y4MWU5LTRjYTMtNjhmOC1iYTQ3LTg3Yzg2NzcxNTQyNQ==; FedAuth=77u/PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0idXRmLTgiPz48U1A+VjE0LDBoLmZ8bWVtYmVyc2hpcHx1cm4lM2FzcG8lM2Fhbm9uIzUyOWY3NjQ1ZWFhOTY1ZjkxNzcyYmJhNWQ3M2E1ZjczMDdkNWY3ZTQzZTgyOTBkODA5NmM2ZmU2NTM5MTY0ZmYsMCMuZnxtZW1iZXJzaGlwfHVybiUzYXNwbyUzYWFub24jNTI5Zjc2NDVlYWE5NjVmOTE3NzJiYmE1ZDczYTVmNzMwN2Q1ZjdlNDNlODI5MGQ4MDk2YzZmZTY1MzkxNjRmZiwxMzM5NDY2MTMwNDAwMDAwMDAsMCwxMzM5NDc0NzQwNTA1ODAyNTksMC4wLjAuMCwyNTgsYjk1ZjUyMjgtNjM5Yy00ZmViLWEwMDQtNzcwYWRjMTQwNjdlLCwsNWM4MjYwMWEtZjI1ZS00YWQxLTg5NDktMmEwMjcyZTIwOWY3LDVjODI2MDFhLWYyNWUtNGFkMS04OTQ5LTJhMDI3MmUyMDlmNyxMWmpoMjVWL1gwV0ZsRUhVSVFuUEpnLDAsMCwwLCwsLDI2NTA0Njc3NDM5OTk5OTk5OTksMCwsLCwsLCwwLCwxOTU4OTcsbWFmVnZaa0RRMTl1WkdIVnhRYWpXNkZLMDBJLCxHb24wS2xHbEEzaVFidjBkcnU4RkJsTU9FdGJWV2tnYW02eW5rZjVQOTkya01pYmYzQk10cTFZcWJNbnN3SExrTVM1bHpIUzhySjM0enlXUG9ydVNYMjJ1TmpjNVoyLytQekZPcWVuUVVnT01YVG1TSUpUVmgvR2dIU2VnKzNlM1lUTjNPWno0Wk1xdDdoM0ZONEJiWTAyMEZQQ3JPSnMvQzlyOFVzMFhZMXVRODlRZUlSTWVjTkI0QWhMMFJlLzNHTS90ZnZDeDV1WnVoWWtiN3lET2g0c2JMRk5uS3FtZC9jWG5jbzFuZi9QcCtuV1FFNUtvRTdTSTc2Sm1VcVcyU1lPWW13TjVpRmlpL093b1hjSklwMWg3bFpUamdVcjFWdmJabzBST0hKMGw1Zmo2bWtXYzdvZWl4RVlNZUJxN09ZUTNQYitVR1pUaUQrY21JbWlKbWc9PTwvU1A+; FeatureOverrides_experiments=[]; MicrosoftApplicationsTelemetryDeviceId=a5ab09b6-1a7f-467e-8b93-10ae64080871; MSFPC=GUID=fc9103de21f54b1fb816e508597e4fa5^&HASH=fc91^&LV=202502^&V=4^&LU=1740503289506; SPA_RT=; ai_session=pWnCfBXn/ClHsgqhSXnoSI^|1750199094486^|1750199094486' \
    -o /gpfs/scratch/nk4167/miniconda/envs/scFoundation_env/models/models1.ckpt