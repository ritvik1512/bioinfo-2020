def matagain (l2,a,b,st):
    print(st)
    ls3 = []
    ls2 = []
    ls3.append(st)
    
    for i in range(0,len(l2)):
        j = l2[i][:]
        ls2.append(j)

    for i in range(0, len(ls2)):
        ls2[i].pop(0)

    for i in range(0,min(len(ls2[a]),len(ls2[b]))):
        distance = round((float(ls2[a][i]) + float(ls2[b][i]))/2,4)
        ls3.append(str(distance))
    ls3.append('0.0')
    l2.append(ls3)
    
    for i in range (0,len(l2)):
        l2[i].append(ls3[i+1])
        if i == len(l2) - 1:
            l2[i].pop(-1)
    
    for i in l2:
        i.pop(a+1)
        i.pop(b+1)
    l2.pop(a)
    l2.pop(b)
    
    return l2 

def minpair (l2):
    while len(l2) > 1:
        a = 0
        b = 0
        min = 1.0
        name = []
        for i in l2:
            name.append(i[0])
            #print(name)
        for i in range(0,len(l2)):
            l3 = []
            l3 = l3 + l2[i]
            l3.pop(0)
            print(l3)     
            for idx in range(0,len(l3)):
                if float(l3[idx]) == 0.0:
                    break
                else:
                    if float(l3[idx]) < min:
                        min = float(l3 [idx])
                        a = i
                        b = idx 
        string = '(' + name[a] +'-'+ name[b] + ')'
        l2 = matagain(l2,a,b,string)

    return string

in_file = open("Distance_data.csv","r")
s1 = in_file.read()

l1 = s1.split(",\n")
#print(l1)
#l1.pop(-1)
l1.pop(0)
l2 = []
name = []

for a in l1:
    i = a[0 : a.find(",")]
    l3 = []
    l3.append(i)
    name = name + l3
    l4 = a[a.find(","):].split(",")
    #print(l4)
    l4.pop(0)
    for a in l4:
        l3.append(a)
    l2.append(l3)

print(len(l2))
print(minpair(l2))
in_file.close()