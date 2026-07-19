import sys,os
from pathlib import Path

script_dir = Path(__file__).resolve().parent
#print(script_dir)
quoted = f'"{script_dir}"'

a = input("Commit message? ")
f = open("gitcommit.bat",mode='w')

f.write('''cd "'''+quoted+'''"
git add .
git commit -m "'''+a+'''"
git push -u origin main''')

f.close()
os.system(quoted+"\\gitcommit.bat")