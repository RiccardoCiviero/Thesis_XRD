def contains_rotate(datfile): # to recognize if I am reading a "normal" omega2theta-scans dat file or a "rotated" omega-scans one 
    with open(datfile) as file:
        for _  in range(3):
            line = file.readline()
            print(line)
        return "rotated" in line
    
datfile = "test.dat"

if contains_rotate(datfile): print("yes")
else: print("no")