import constants
import functools
import collections
import logging
import numpy as np

class Cacher(object):
    """
    кешируем в sqlite то, что можем
    """
    def __init__(self):
        if constants.sql_enabled: # если кеширование в БД есть
            import sqlite3
            self.db = sqlite3.connect('cache.db')
            self.cur = self.db.cursor()
            self.cur.execute("CREATE TABLE IF NOT EXISTS TimeFromTempreture(value REAL NULL, tempeture REAL);")
            self.cur.execute("CREATE TABLE IF NOT EXISTS TempretureFromTime(value REAL NULL, tempeture REAL);")
            self.cur.execute("CREATE TABLE IF NOT EXISTS DerriviateTdtFromTempreture(value REAL NULL, tempeture REAL);")
            self.cur.execute("CREATE TABLE IF NOT EXISTS TempreNuFromTempreture(value REAL NULL, tempeture REAL);")
            self.cur.execute("CREATE TABLE IF NOT EXISTS LambdaNPFromTempr(value REAL NULL, tempeture REAL);")
            self.cur.execute("CREATE TABLE IF NOT EXISTS LambdaPNFromTempr(value REAL NULL, tempeture REAL);")

        self._INNER_CACHE_SIZE_ = 100
        self._inner_cache_ = collections.defaultdict(dict)

    def sql_base_cache(self, func):
        """
        кеширование температуры фотонов, нейтрино и времени
        """
        @functools.wraps(func)
        def inner(*args, **kwargs):
            tab_name = constants.func_name_to_db_name[func.__name__]
            coef = 1.0
            cached = self.cur.execute("SELECT * FROM {table} WHERE tempeture={arg};".format(
                arg=args[0]*coef, table=tab_name)
                )
            res = self.cur.fetchone() # ищем результат в кеше
            if res: # если нашли температуру
                return res[0]
            if not res: # не нашли - вставим
                res = func(*args, **kwargs)
                self._inner_cache_[tab_name][args[0]*coef] = res
                if len(self._inner_cache_[tab_name]) > self._INNER_CACHE_SIZE_:
                    for tem, value in self._inner_cache_[tab_name].items():
                        exe_str = "INSERT OR REPLACE INTO {table} (value, tempeture) VALUES ({valu}, {tem});".format(
                            table=tab_name, tem=tem, valu=value
                            )
                        self.cur.execute(exe_str)
                    self._inner_cache_[tab_name] = dict()
                    import sqlite3
                    try:
                        self.db.commit()
                    except sqlite3.OperationalError as e:
                        logging.error("db is lock")
                return res

        return inner if constants.sql_enabled else func

    def __del__(self):
        if constants.sql_enabled:
            for tab_name, di in self._inner_cache_.items():
                for tem, value in di.items():
                    exe_str = "INSERT OR REPLACE INTO {table} (value, tempeture) VALUES ({valu}, {tem});".format(
                        table=tab_name, tem=tem, valu=value
                        )
                    self.cur.execute(exe_str)
            self.db.commit()
            self.db.close()


cacher = Cacher()

