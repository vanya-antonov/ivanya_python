#!/usr/bin/env python3

import sys, os, re, shutil
import MySQLdb
import MySQLdb.cursors
from pprint import pprint
from random import randint

from mylib.local import DB_ACCOUNTS

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
	# https://dev.mysql.com/doc/connector-python/en/connector-python-api-mysqlcursor-execute.html
	def exec_sql_ar(self, sql, *params):
		self.cur.execute(sql, params)
		return self.cur.fetchall()
	
	#---
	def exec_sql_nr(self, sql, *params):
		self.cur.execute(sql, params)
	
	#---
	# Useful for INSERT statements: the sql should be truncated like 'INSERT INTO users (id, name, descr)'
	def exec_sql_in(self, sql, *params):
		# Make a list by repearing the '%s' string: https://stackoverflow.com/a/3459131/310453
		tail_l = ['%s'] * len(params)
		sql += ' VALUES (' + ','.join(tail_l) + ')'
		self.cur.execute(sql, params)
	
	#---
	def get_random_db_id(self, table, field):
		sql = 'SELECT 1 AS count FROM ' + table + ' WHERE ' + field + '=%s'
		num = randint(10**8, 10**9 - 1)
		while len(self.exec_sql_ar(sql, num)) > 0:
			num = randint(10**8, 10**9 - 1)
		return num

if __name__ == '__main__':
	print('Hello!!!')
	budb = BaseUtilDB(db_account='gtdb2')
	all_rows = budb.exec_sql_ar("SELECT * FROM users WHERE id = %s", -2)
	pprint(all_rows)
	


