	


class colors:
	"""
    reset all colors with colors.reset
    two subclasses fg for foreground and bg for background.
    use as colors.subclass.colorname.
    i.e. colors.fg.red or colors.bg.green
    also, the generic bold, disable, underline, reverse, strikethrough,
    and invisible work with the main class
    i.e. colors.bold
	"""
	reset='\033[0m'
	bold='\033[01m'
	disable='\033[02m'
	underline='\033[04m'
	reverse='\033[07m'
	strikethrough='\033[09m'
	invisible='\033[08m'
	class fg:
		black='\033[30m'
		red='\033[31m'
		green='\033[32m'
		orange='\033[33m'
		blue='\033[34m'
		purple='\033[35m'
		cyan='\033[36m'
		lightgrey='\033[37m'
		darkgrey='\033[90m'
		lightred='\033[91m'
		lightgreen='\033[92m'
		yellow='\033[93m'
		lightblue='\033[94m'
		pink='\033[95m'
		lightcyan='\033[96m'
	class bg:
		black='\033[40m'
		red='\033[41m'
		green='\033[42m'
		orange='\033[43m'
		blue='\033[44m'
		purple='\033[45m'
		cyan='\033[46m'
		lightgrey='\033[47m'
		
		
		

class Screen:
	"""
	Handles the output screen in NUDGE
	
	The top header is always displayed
	
	"""
	
	def __init__(self):
		self.lines = 24
		self.bottom_buffer = 1
		self.width = 80
		self.default_top = 3
		self.current_line = 3
		self.header = '----------------- NUDGE: NUclear Database GEneration software -----------------'
		return
	
	def PaintBackground(self, color, lines=None, width=None):
		if lines == None: lines = self.lines
		if width == None: width = self.width
		for i in range(lines):
			print(color,' '*width)

	def PrintAt(self, text, x=1, y=None, color=colors.fg.green):
		if y == None: y = self.current_line
		x = int(x)
		y = int(y)
		if x >= self.width: x = 255
		if y >= self.lines: y = 255
		if x <= 0: x = 0
		if y <= 0: y = 0
		
		HORIZ = str(x)
		VERT = str(y)
		
		# Plot the text at the starting at position HORIZ, VERT...
		print(color + "\033["+VERT+";"+HORIZ+"f" + text)
		self.current_line += 1
	
	# Initializes the screen (with correct background
	def InitScreen(self):
		self.PaintBackground(colors.bg.black) # Paints the background color
		self.Header()
	
	# Prints the header
	def Header(self):
		self.PrintAt(self.header, 0,0, colors.fg.cyan)
		print(colors.fg.green)
		
	# Prints the help screen
	def HelpScreen(self):
		help_text = 'NUDGE is used to generate xsgen libraries using scout runs and global surrogate modeling' \
					+ '\n -g (guided): start NUDGE in guided mode \n -d (database): used for the database path \n -h (help):  help screen'
		print(help_text)
	
	def UpdateInfo(self, database):
		self.complete_slibs = str(database.complete_slibs)
		self.slibs = str(len(database.slibs))
		self.complete_flibs =  str(database.complete_flibs)
		self.flibs = str(len(database.flibs))
		text = 'Scout libs: ' + self.complete_slibs + '/' + self.slibs + \
				'  -  Full libs: ' + self.complete_flibs + '/' + self.flibs + \
				'  -  Dimensions: ' + str(database.dimensions)
		self.PrintAt(text, 2, 2, color=colors.fg.lightblue)
		text = 'Varied inputs: ' + str(database.varied_ips)
		self.PrintAt(text, 2, 3, color=colors.fg.lightblue)
	
	def Clean(self):
		self.InitScreen()
		self.UpdateInfo()
		
		
		
		
		
		
		
		
		
		
		
		
