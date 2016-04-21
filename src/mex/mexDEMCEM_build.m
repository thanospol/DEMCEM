% build DIRECT_WS
mex -output mexDIRECT_WS_ST_RWG    -O DIRECT_WS_ST_RWG_mex.cpp   DIRECT_WS_ST_RWG.cpp 
mex -output mexDIRECT_WS_EA_RWG    -O DIRECT_WS_EA_RWG_mex.cpp   DIRECT_WS_EA_RWG.cpp
mex -output mexDIRECT_WS_VA_RWG    -O DIRECT_WS_VA_RWG_mex.cpp   DIRECT_WS_VA_RWG.cpp 
% build DIRECT_SS
mex -output mexDIRECT_SS_EA_RWG    -O DIRECT_SS_EA_RWG_mex.cpp   DIRECT_SS_EA_RWG.cpp 
mex -output mexDIRECT_SS_VA_RWG    -O DIRECT_SS_VA_RWG_mex.cpp   DIRECT_SS_VA_RWG.cpp 
mex -output mexDIRECT_SS_EA_nxRWG  -O DIRECT_SS_EA_nxRWG_mex.cpp DIRECT_SS_EA_nxRWG.cpp 
mex -output mexDIRECT_SS_VA_nxRWG  -O DIRECT_SS_VA_nxRWG_mex.cpp DIRECT_SS_VA_nxRWG.cpp 