import pacal

def sim(g, up, dw):
	x = pacal.BinomialDistr(g, up)
	y = pacal.BinomialDistr(g, dw)
	return x-y
