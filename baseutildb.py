#!/usr/bin/env python3

import MySQLdb
import MySQLdb.cursors
from local import DB_ACCOUNTS
from pprint import pprint

class BaseUtilDB:
	def __init__(self, db_account="default"):
		db_info = DB_ACCOUNTS.get(db_account)
		if db_info is None:
			raise Exception("Unknown db_account = '%s'" % db_account)
		
		self.db = MySQLdb.connect(host   = db_info['host'],  # your host 
		                          db     = db_info['name'],  # name of the database
		                          user   = db_info['user'],  # username
		                          passwd = db_info['pass'],  # password
		)
		
		# https://stackoverflow.com/a/2180257/310453
		self.cur = self.db.cursor(MySQLdb.cursors.DictCursor)
	
	#---
	def exec_SQL_ar(self,sql):
		self.cur.execute(sql)
		return self.cur.fetchall()

if __name__ == '__main__':
	print('Hello!!!')
	budb = BaseUtilDB(db_account='ivdb')
	all_rows = budb.exec_SQL_ar("SELECT * FROM users limit 2")
	pprint(all_rows)
	


