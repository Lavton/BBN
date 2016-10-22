import constants
import functools
# sql_enabled, func_name_to_db_name
class Cacher(object):
    """docstring for Cacher"""
    def __init__(self):
        if constants.sql_enabled: # если кеширование в БД есть
            import sqlite3
            import collections
            self.db = sqlite3.connect('cache.db')
            self.cur = self.db.cursor()
            self.cur.execute("CREATE TABLE IF NOT EXISTS TimeFromTempreture(value REAL NULL, tempeture REAL);")
            self.cur.execute("CREATE TABLE IF NOT EXISTS TempreNuFromTempreture(value REAL NULL, tempeture REAL);")

            self._INNER_CACHE_SIZE_ = 100
            self._inner_cache_ = collections.defaultdict(dict)
            # self._num_of_exect
        # if constants.sql_enabled:
            # print("initialize")
        # super(Cacher, self).__init__()

    def sql_tempreture_cache(self, func):
        """
        кеширование температуры фотонов, нейтрино и времени
        """
        @functools.wraps(func)
        def inner(*args, **kwargs):
            tab_name = constants.func_name_to_db_name[func.__name__]
            coef = 1.0
            if kwargs.get("units", "") == "K":
                coef = constants.k_b # если температура была дана в Кельвинах, переводин в eV
            res = self.cur.fetchone() # ищем результат в кеше
            if res: # если нашли температуру
                return res[0]
            if not res: # не нашли - вставим
                res = func(*args, **kwargs)
                exe_str = "INSERT INTO {table} (value, tempeture) VALUES ({valu}, {tem});".format(
                    table=tab_name, tem=args[0]*coef, valu=res
                )
                self.cur.execute(exe_str)
                self.db.commit()
                return res

        return inner if constants.sql_enabled else func


    def __del__(self):
        if constants.sql_enabled:
            self.db.close()


# cacher = Cacher()
# sql_tempreture_cache = cacher.sql_tempreture_cache
# import sqlite3
# import collections

# db = sqlite3.connect('cache.db')
# cur = db.cursor()

