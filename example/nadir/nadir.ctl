# ======================================================================
# Forward model...
# ======================================================================

# Table directory...
TBLBASE = ./airs

# Emitters...
NG = 1
EMITTER[0] = CO2

# Channels...
ND = 3
NU[0] = 667.7820
NU[1] = 668.5410
NU[2] = 669.8110

# Output...
WRITE_BBT = 1

# use the GPU: 0:never, 1:always, -1:if possible
USEGPU = -1
WRITE_BINARY = 0
READ_BINARY = 0
