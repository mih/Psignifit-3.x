#
# Regular cron jobs for the psignifit3 package
#
0 4	* * *	root	[ -x /usr/bin/psignifit3_maintenance ] && /usr/bin/psignifit3_maintenance
