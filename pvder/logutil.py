import os
import json

from pvder.defaults import logConfig
from pvder.exceptionutil import ExceptionUtil


LogUtil=ExceptionUtil()
baseDir=os.path.dirname(os.path.abspath(__file__))
logConfig['logFilePath']=os.path.join(baseDir,logConfig['logFilePath'])
LogUtil.create_logger(**logConfig)
