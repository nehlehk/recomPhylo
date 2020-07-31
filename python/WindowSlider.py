
a =[1,2,1,2,10,2,1,2,1,2,3,4,5,6,2]

def slidingWindow(sequence,winSize,step=1):
    # Verify the inputs
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")

    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence)-winSize)/step)+1

    # Do the work
    for i in range(0,numOfChunks*step,step):
        yield sequence[i:i+winSize]


result = slidingWindow(a,3,1)

print(result)
