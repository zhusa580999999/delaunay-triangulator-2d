ADD_LIBRARY(reader SHARED
            spreader_spa.cpp
	    spreader_spb.cpp
	    spreader_node.cpp
	    spreader_raw.cpp
	    spreader_raw_d.cpp
	    spreader_tiles.cpp
	    smreader_sma.cpp
	    smreader_smb.cpp
	    svreader_sva.cpp
	    svreader_svb.cpp)

ADD_LIBRARY(writer SHARED
            spwriter_spa.cpp
	    spwriter_spb.cpp
	    spwriter_raw.cpp
	    smwriter_sma.cpp
	    smwriter_smb.cpp)

# C interface to reader and writer
ADD_LIBRARY(cio SHARED cio.cpp)

TARGET_LINK_LIBRARIES(cio reader writer)
