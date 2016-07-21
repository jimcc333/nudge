	


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
		
		
		

class screen:
	"""
	Handles the output screen in NUDGE
	
	The top header is always displayed
	
	"""
	lines = 24
	width = 80
	header = '----------------- NUDGE: NUclear Database GEneration software -----------------'
	
	def PaintBackground(color, lines=lines, width=width):
		for i in range(lines):
			print(color,' '*width)

	def PrintAt(text, x=1, y=5, color=colors.fg.green):
		x = int(x)
		y = int(y)
		if x >= 255: x = 255
		if y >= 255: y = 255
		if x <= 0: x = 0
		if y <= 0: y = 0
		
		HORIZ = str(x)
		VERT = str(y)
		
		# Plot the text at the starting at position HORIZ, VERT...
		print(color + "\033["+VERT+";"+HORIZ+"f" + text)	
	
	def InitScreen():
		screen.PaintBackground(colors.bg.black) # Paints the background color
		screen.PrintAt(screen.header, 0,0, colors.fg.cyan) # Places the header
		print(colors.fg.green)
		
	def HelpScreen():
		help_text = 'NUDGE is used to generate xsgen libraries using scout runs and global surrogate modeling' \
					+ '\n -g (guided): start NUDGE in guided mode \n -d (database): used for the database path \n -h (help):  help screen'
		print(help_text)
	
	def UpdateInfo(database):
		text = 'Scout libs: ' + str(database.complete_slibs) + '/' + str(len(database.slibs)) + \
				'  -  Full libs: ' + str(database.complete_flibs) + '/' + str(len(database.flibs)) + \
				'  -  Dimensions: ' + str(database.dimensions)
		screen.PrintAt(text, 2, 2)
		
		
		
		
		
		
		
		
		
		
		
		
