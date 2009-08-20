include $(top_srcdir)/config/MakeGlobal.am

lib_LIBRARIES = libMADmisc.a

bin_PROGRAMS = interp3

thisincludedir = $(includedir)/misc/
thisinclude_HEADERS = misc.h ran.h \
                      muParserBase.h muParserBytecode.h muParserCallback.h \
                      muParserDef.h muParserDLL.h muParserError.h \
                      muParserFixes.h muParser.h muParserInt.h \
                      muParserStack.h muParserTest.h muParserToken.h \
                      muParserTokenReader.h phandler.h


libMADmisc_a_SOURCES = checksum_file.cc position_stream.cc gprofexit.cc ran.cc \
                        \
                       misc.h ran.h \
                      muParserBase.cpp muParserBytecode.cpp \
                      muParserCallback.cpp \
                      muParser.cpp muParserError.cpp muParserTokenReader.cpp \
                      muParserBase.h muParserBytecode.h muParserCallback.h \
                      muParserDef.h muParserDLL.h muParserError.h \
                      muParserFixes.h muParser.h muParserInt.h \
                      muParserStack.h muParserTest.h muParserToken.h \
                      muParserTokenReader.h phandler.h

interp3_SOURCES = interp3.cc
