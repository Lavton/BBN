Код моделирует образование элементов во время первичного нуклеосинтеза.
Запуск - `python3 nucleosynthesis.py`


Логика: все слогаемые с реакцией между элементами записаны в файле с самым тяжёлым элементом в elements/*

Фунции объединяются в одну (см. https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html ) в elements/register. И уже технические изменения этих функций, ровно как и сам подсчёт, происходит в nucleosynthesis