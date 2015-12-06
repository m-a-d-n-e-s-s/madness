include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

lib_LTLIBRARIES = libMADmuparser.la

thisincludedir = $(includedir)/madness/external/muParser
thisinclude_HEADERS = muParserBase.h muParserBytecode.h muParserCallback.h \
                      muParserDef.h muParserDLL.h muParserError.h \
                      muParserFixes.h muParser.h muParserInt.h \
                      muParserStack.h muParserTest.h muParserToken.h \
                      muParserTokenReader.h


libMADmuparser_la_SOURCES = muParserBase.cpp muParserBytecode.cpp \
                      muParserCallback.cpp \
                      muParser.cpp muParserError.cpp muParserTokenReader.cpp \
                      muParserBase.h muParserBytecode.h muParserCallback.h \
                      muParserDef.h muParserDLL.h muParserError.h \
                      muParserFixes.h muParser.h muParserInt.h \
                      muParserStack.h muParserTest.h muParserToken.h \
                      muParserTokenReader.h
libMADmuparser_la_LDFLAGS = -version-info 0:0:0

libraries:	$(lib_LTLIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
