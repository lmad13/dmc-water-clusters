def printArray(a):
    string=[]
    for row in a:
        string.append(str(row).replace('[',' ').replace(']',' '))
    return string
    
