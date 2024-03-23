def showAttributes(att,filter="none"):
    for element in att:
        if filter=="none":
            print(element)
        elif element.startswith(filter):
            print(element)

def getAttributes(obj,filter="none",showHidden=0):
    attributes = []
    for d in dir(obj):
        if showHidden==0: # dont show hidden variables
            if not d.startswith("_") and not(d == "custom_attribute"):
                attributes.append(d)
        elif showHidden==1: # show level 1 hidden variables
            if not d.startswith("__") and not(d == "custom_attribute"):
                attributes.append(d)
        else:
            attributes.append(d)
    
    showAttributes(attributes,filter)
    
    return attributes