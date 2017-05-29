import constants
import functools
import collections
import logging
import numpy as np

# sql_enabled, func_name_to_db_name
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
            # if kwargs.get("units", "") == "K":
                # coef = constants.k_b # если температура была дана в Кельвинах, переводин в eV
            cached = self.cur.execute("SELECT * FROM {table} WHERE tempeture={arg};".format(
                arg=args[0]*coef, table=tab_name)
                )
            res = self.cur.fetchone() # ищем результат в кеше
            if res: # если нашли температуру
                # logging.debug("cache from SQL")
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
                    try:
                        self.db.commit()
                    except sqlite3.OperationalError as e:
                        logging.error("db is lock")
                return res

        return inner if constants.sql_enabled else func

    # def np_array_to_list_decor(self, func):
    #     """
    # Кеширование
    #     """

    #     @functools.wraps(func)
    #     def inner(*args, **kwargs):
    #         new_args = []
    #         new_kwargs = dict()
    #         for i in range(len(args)):
    #             if isinstance(args[i], np.ndarray):
    #                 print(tuple(args[i]))
    #                 new_args.append(tuple(args[i]))
    #             else:
    #                 new_args.append(args[i])
    #         for k, v in kwargs.items():
    #             if isinstance(v, np.ndarray):
    #                 new_kwargs[k] = tuple(v)
    #             else:
    #                 new_kwargs[k] = v
    #         return functools.lru_cache(func)(*args, **kwargs)
    #     return inner


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

