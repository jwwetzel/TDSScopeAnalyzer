CC			=	g++
CFLAGS		=	-c -g -Wall `root-config --cflags`
CXX			=	`root-config --cxx`
CXXFLAGS	=	`root-config --cflags`
LDFLAGS		=	`root-config --ldflags`
LDLIBS		=	`root-config --glibs`

SOURCES		=	TDSwfmAnalyzer.cxx
HEADERS		=	TDSwfmAnalyzer.h LinkDef.h
OBJECTS		=	$(SOURCES:.cxx=.o)
EXECUTABLE	=	TDSwfmAnalyzer
TDSwfmAnalyzer : TDSwfmAnalyzer.cxx

$(EXECUTABLE): $(SOURCES)
	$(CXX) $(CXXFLAGS) -W -Wall -o $@ $^ $(LDLIBS)

clean:
	rm ./*~ ./*.o
