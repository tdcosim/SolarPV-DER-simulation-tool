import os
import sys
import logging
import linecache
import pdb
import json
from dateutil import parser


#===================================================================================================
#===================================================================================================
#===================================================================================================
class ExceptionUtil(object):
	"""
	Sample usage:
	=============
	from exceptionutil import ExceptionUtil

	class Foo(ExceptionUtil):
		def __init__(self):
			self.create_logger("sampleLogger",logFilePath='path/to/myLog.log',logLevel=logging.DEBUG,
			mode='w')
			return None

		def bar(self,a,b):
			try:
				code
				self.logger.log(level=10,msg='')
				some more code
			except:
				self.exception_handler()
	"""

#===================================================================================================
	def __init__(self):
		super(ExceptionUtil,self).__init__()
		self._defaultFormatterStr='%(asctime)s::%(name)s::%(filename)s::%(funcName)s::'+\
		'%(levelname)s::%(message)s::%(threadName)s::%(process)d'
		self._defaultFormatterStrSep='::'
		self._formatterStr=None
		self._formatterStrSep=None
		return None

#===================================================================================================
	def create_logger(self,loggerName,logFilePath,logLevel,formatterStr=None,formatterStrSep=None,mode='a'):
		try:
			# create logger
			self.logger = logging.getLogger(loggerName)
			self._logFilePath=logFilePath
			
			# create formatter and add it to the handlers
			if not formatterStr:
				formatterStr=self._defaultFormatterStr
				formatterStrSep=self._defaultFormatterStrSep

			self._formatterStr=formatterStr
			self._formatterStrSep=formatterStrSep
			self.set_logger(logLevel,mode)

		except:
			raise

#===================================================================================================
	def set_logger(self,logLevel,mode):
		try:
			self.logger.setLevel(logLevel)
			# create file handler which logs messages
			fh = logging.FileHandler(self._logFilePath,mode=mode)

			formatter = logging.Formatter(self._formatterStr)
			fh.setFormatter(formatter)
			fh.setLevel(logLevel)

			self.logger.addHandler(fh)
		except:
			raise

#===================================================================================================
	def _get_exception(self):
		try:
			res={}
			exc_type, exc_obj, tb = sys.exc_info()
			f = tb.tb_frame
			res['lineno'] = tb.tb_lineno
			res['filename'] = f.f_code.co_filename
			res['funcName']=f.f_code.co_name
			linecache.checkcache(res['filename'])
			res['line'] = linecache.getline(res['filename'], res['lineno'], f.f_globals).strip()
			res['reason']='{}'.format(exc_obj)
			return res
		except:
			raise

#===================================================================================================
	def exception_handler(self,additionalInfo=None):
		try:
			err=self._get_exception()
			if additionalInfo:
				err['additionalInfo']=additionalInfo
			self.logger.error(json.dumps(err))
			raise
		except:
			raise

#===================================================================================================
	def log(self,level=None,msg=None):
		try:
			if not level:
				level=logging.ERROR
			if level==logging.ERROR:# will also raise
				self.exception_handler(msg)
			else:
				self.logger.log(level,msg)
		except:
			raise

#===================================================================================================
	def get_logs(self,filterResult=None,logFilePath=None,formatterStr=None,formatterStrSep=None):
		"""
		Sample call:
		============
		self.get_logs(logFilePath='foo.log')
		"""
		try:
			if not formatterStr and self._formatterStr:
				formatterStr=self._formatterStr
			elif not formatterStr and not self._formatterStr:
				formatterStr=self._defaultFormatterStr

			if not formatterStrSep and self._formatterStrSep:
				formatterStrSep=self._formatterStrSep
			elif not formatterStrSep and not self._formatterStrSep:
				formatterStrSep=self._defaultFormatterStrSep

			if not logFilePath:
				logFilePath=self._logFilePath
			f=open(logFilePath);
			logs=f.read(); f.close()

			res=[]; addThisItem=True
			header=[entry for entry in formatterStr.split(formatterStrSep)]
			for thisLine in logs.splitlines():
				thisLog={}
				for item,val in zip(header,thisLine.split(formatterStrSep)):
					thisLog[item.strip()]=val.strip()
				if '%(asctime)s' in thisLog:# convert to datetime object
					thisLog['%(asctime)s']=parser.parse(thisLog['%(asctime)s'])
				if '%(levelname)s' in thisLog and thisLog['%(levelname)s']=='ERROR':
					thisLog['%(message)s']=json.loads(thisLog['%(message)s'])
					if '%(filename)s' in thisLog:
						thisLog['%(filename)s']=os.path.basename(thisLog['%(message)s']['filename'])
					if '%(funcName)s' in thisLog:
						thisLog['%(funcName)s']=os.path.basename(thisLog['%(message)s']['funcName'])

				if filterResult:
					if isinstance(filterResult,dict):
						addThisItem=True
						for entry in filterResult:
							if not entry in thisLog or thisLog[entry]!=filterResult[entry]:
								addThisItem=False
								break
				if addThisItem:
					res.append(thisLog)

			return res
		except:
			raise


