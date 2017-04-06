#Erik Sandberg
#Look Say Sequence

def where_not_eq(string_array,value):
    for i, val in enumerate(string_array):
        if int(val) != int(value):
            return i

array = ['1']
steps = 7
for i in range(steps):
    string = array[-1]
    value = string[0]
    unfinished = ''
    counter = 0
    while True:
        counter +=1 
        summed = where_not_eq(string,value)
        if summed == None:
            summed = 1 + len(string[int(value):])
        unfinished += str(summed) + value
        try:
            value = string[summed]
        except:
            pass
        string = string[summed:]
        if len(string) == 0:
            break
    array.append(unfinished)

print array
