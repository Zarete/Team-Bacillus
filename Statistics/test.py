# -*- coding: utf-8 -*-

import Stats as stats

seq = stats.getRandomSequence(1000)

test = 'hello'

print(test[0])

print(seq)
print(len(seq['data']))

print(stats.shuffle(seq))

stats.writeCSV('test', ',', seq)