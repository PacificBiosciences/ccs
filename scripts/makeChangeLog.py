""" File to generate a CHANGELOG from git commits.  This will parse the git log
command and convert any commit with `**` in the title into a changelog event """

cl = open("CHANGELOG", 'w')
cl.write("Recent Program Changes:\n\n")

import os
from subprocess import *

GIT_COMMIT_FIELDS = ['id', 'author_name', 'author_email', 'date', 'message']
GIT_LOG_FORMAT = ['%H', '%an', '%ae', '%ad', '%s']

p = Popen('git log --format="%s"' % GIT_LOG_FORMAT, shell=True, stdout=PIPE)
(log, _) = p.communicate()
log = log.strip().split("\n")

log = [row.strip().replace("'","").split(",") for row in log]
log = [dict(zip(GIT_COMMIT_FIELDS, row)) for row in log]

for commit in log:
	msg = commit["message"]
	if msg.count("**") > 0:
		dt = commit["date"].split(" ")
		date = "-".join(dt[2:4])
		msg = msg.replace("]","").replace("**","")
		id = commit["id"][1:7]
		cl.write("\t*" + " " + msg + " " + date + " " + id + "\n")

cl.write("\n\nFile is automatically created by makeChangeLog.py\n\n\n")

cl.close()