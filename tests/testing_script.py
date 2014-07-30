from __future__ import print_function
import subprocess
import os
import sys
import time



UNKNOWN = 0
WHITE = 1
BLACK = 2

HEX = 1
REX = 2
CHEX = 3
Y = 4
GALE = 5



#Send cmd to the terminal for it to execute
def external_command(cmd): 
    process = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE)

    # wait for the process to terminate
    out, err = process.communicate()
    errcode = process.returncode

    return errcode, out, err
    
    
class Test:
	
	testCount = 0
   

	def __init__(self, test_number, expected_outcome, passed_MCTS, passed_PNS):
		self.test_number= test_number
		self.expected_outcome = expected_outcome
		self.passed_MCTS = passed_MCTS
		self.passed_PNS = passed_PNS

		Test.testCount += 1
   
	def displayTestCount(self):
		print ("Total Tests %d" % Test.testCount)
		
def getTests(tests):
	f = open(tests)
	t = f.readlines()
	t = [x.strip('\n') for x in t]
	return t

#Reads from the test file created by user and parses out the relevant data
def getCommands(tests, game_type):
	#This is a check to see if we have found a test
	found_test = False
	#Check to see if we have found the commands for that test
	found_commands = False
	found_size = False

	#List of Commands to execute in Rex Program
	commands = []
	test_commands = []
	
	theoretical_results = []

	test_number = UNKNOWN
	expected_outcome = UNKNOWN
	#Iterate through the test file and
	#parse out the commands
	for r in range(0, len(tests)):
		if ('test' == tests[r]):
			found_test = True
			test_number = tests[r+1]
			continue
		elif ('expected_outcome' == tests[r]):
			expected_outcome = tests[r+1]
			continue
		elif ('commands' == tests[r]):
			found_commands = True
			continue
		#Else if we have come to the end of a test case
		elif ('endtest' == tests[r]):
			if (found_size == False):
				print("Error found in test script.  A test is missing the 'size' command")
			found_size = False
			found_test = False
			found_commands = False
			#Append either mcts or pns or solver to the commands list
			test_commands.append("mcts")
			test_commands.extend(commands)
			if (game_type == GALE):
				test_commands.append("solver")
				test_commands.extend(commands)
			else:
				test_commands.append("pns")
				test_commands.extend(commands)
			commands = []
			#Append a test to the list with the specified test number and
			#expected outcome.  The last two are False because we do not 
			#know if the test passed yet
			theoretical_results.append(Test(str(test_number) + '_mcts', expected_outcome, False, False))
			if (game_type == GALE):
				theoretical_results.append(Test(str(test_number) + '_solver', expected_outcome, False, False))
			else :
				theoretical_results.append(Test(str(test_number) + '_pns', expected_outcome, False, False))
			continue
		#Append all of the commands to the list
		if (found_commands):
			commands.append(tests[r])
		if ('size' in tests[r]):
			found_size = True
	
	test_commands.append("quit")			
	return test_commands, theoretical_results

#Analyze the actual results of the script versus what should have happened
def analyzeResults(actual_results, theoretical_results, game_type):
	t = 0		
	
	for a in actual_results:
		print (a)
			
	if(len(actual_results) != len(theoretical_results)):
		print("Error in analysis: Actual Results #: "  + str(len(actual_results)) + "  Theoretical Results #: " + str(len(theoretical_results)))
	
		
	number_of_tests = len(theoretical_results)
	number_of_passed_tests = 0
	number_of_MCTS_passed_tests = 0
	number_of_PNS_passed_tests = 0
	number_of_solver_passed_tests = 0
	number_of_agreeing_tests = 0
	#List of all failed tests
	failed_tests = []
	#List of all conflicting tests between MCTS and PNS or Solver
	conflicting_tests = []
	unknown_tests = []
	for a in range(0, len(actual_results)):
		
		#Check that the actual result matches the theoretical result
		if (("Solved as a win" in actual_results[a] and (theoretical_results[t].expected_outcome is str(1))) \
		or ("Solved as a loss" in actual_results[a] and (theoretical_results[t].expected_outcome is str(0)))):
			number_of_passed_tests += 1
			if('_mcts' in theoretical_results[t].test_number):
				number_of_MCTS_passed_tests += 1
			elif ('_pns' in theoretical_results[t].test_number):
				number_of_PNS_passed_tests += 1
			elif ('_solver' in theoretical_results[t].test_number):
				number_of_solver_passed_tests += 1
		#If they did not match it is a failed test
		elif("Unknown" in actual_results[a]):
			unknown_tests.append(theoretical_results[t])
		else:
			failed_tests.append(theoretical_results[t])
		
		#If this was an MCTS test, then the next one is an identical test but for PNS
		#Check if that one has the same result
		if('_mcts' in theoretical_results[t].test_number):
			if(actual_results[a] == actual_results[a + 1]):
				number_of_agreeing_tests += 2
			#Else if they did not agree, note it
			else:
				conflicting_tests.append(theoretical_results[t])
				conflicting_tests.append(theoretical_results[t+1])
		t += 1
			
	print("\n\n*********Analysis***********\n")
	print("Number of Tests: " + str(number_of_tests))
	print(str(number_of_passed_tests) + "/" + str(number_of_tests) + " Tests Passed")
	print(str(number_of_agreeing_tests) + "/" + str(number_of_tests) + " Tests agree with each other")
	print(str(number_of_MCTS_passed_tests) + "/" + str(number_of_tests / 2) + " MCTS Tests Passed")
	if (game_type == GALE):
		print(str(number_of_solver_passed_tests) + "/" + str(number_of_tests / 2) + " Solver Tests Passed")
	else:
		print(str(number_of_PNS_passed_tests) + "/" + str(number_of_tests / 2) + " PNS Tests Passed")
	print(str(len(failed_tests)) + "/" + str(number_of_tests) + " Tests Failed")
	print(str(len(unknown_tests)) + "/" + str(number_of_tests) + " Tests were unknown")
	print("Failed Tests: " )
	for f in failed_tests:
		print (f.test_number)
	print("Conflicting Tests: ")
	for c in conflicting_tests:
		print(c.test_number)
	print("Unknown Tests: ")
	for u in unknown_tests:
		print(u.test_number)
		
def parseResults(results):
	game_results = open(results, 'r')
	for r in game_results:
		actual_results = r.split('\\n')
	#Check to see if everything found a solution
	#If not, tag an 'Unknown' onto it
	found_solution = False
	for r in range(0, len(actual_results)):
		if('Solved as a' in actual_results[r]):
			found_solution = True
		if ('PV:' in actual_results[r]):
			if (not found_solution):
				actual_results[r] = "Unknown"
			found_solution = False	
		if ("It is the other player's turn" in actual_results[r] or "Unknown command" in actual_results[r]):
			print("Error found on line " + str(r))
			
	
	#Remove all lines from the results data except for the lines
	#that contain relevant information (i.e. Whether it solved it as a win
	#or a loss or whether it was unknown)				
	actual_results[:] = [a for a in actual_results if not ("Solved as a" not in a) or not ("Unknown" not in a)]
	return actual_results
	
def runTests(tests, game_type):
	if (game_type == REX):
		game = "Rex"
	elif (game_type == HEX):
		game = "Hex"
	elif (game_type == CHEX):
		game = "Cylindrical Hex"
	elif (game_type == Y):
		game = "Y"
	elif (game_type == GALE):
		game = "Bridg It"
	
	#Parse out the commands and the results from the test file
	#to obtain a list of commands that the program can execute
	test_commands, theoretical_results = getCommands(tests, game_type)
	
	#Print the commands to a file
	#which will be used to later to
	#pipe into the game
	cf = open('commands_file', 'w')
	for t in test_commands:	
		print(t, file=cf)
	cf.close()	
	
	#Clear the terminal to make it look nicer
	os.system('clear')
	print("Performing " + game + " Tests")
	
	#Redirect Stdout to a file temporarily 
	temp = sys.stdout
	#Create a file and write the test results to it
	sys.stdout = open('test_results', 'w')
	
	#Start a timer
	start = time.clock()
	
	#Run the program and pipe in the commands
	#Save the output to a file
	if (game_type == REX):
		output = external_command("../trex -f commands_file")
	elif (game_type == HEX):
		output = external_command("../nhex -f commands_file")
	elif (game_type == CHEX):
		output = external_command("../chex -f commands_file")
	elif (game_type == Y):
		output = external_command("../moy -f commands_file")
	elif (game_type == GALE):
		output = external_command("../gale -f commands_file")
	print(output)
	#Stop the timer
	elapsed = (time.clock() - start)
	#Set stdout back to what it was normally
	sys.stdout = temp
	
	print("Completed " + game + " Tests.  Time Taken: " + str(elapsed) + "s. \n\n Analyzing Results...")
	#Parse Results
	actual_results = parseResults('test_results')
	
	#Analyze the results
	analyzeResults(actual_results, theoretical_results, game_type)

	#Delete unnecessary files
	#os.remove("test_results")
	#os.remove("commands_file")

	
	

if __name__ == '__main__':
	#Open Test File
	game = raw_input("Which program do you want to run tests for? (nhex, rex, chex, gale): ")
	if ("rex" in game):
		rex_tests = getTests('rex_tests')
		runTests(rex_tests, REX)
	elif("nhex" in game):
		hex_tests = getTests('hex_tests')
		runTests(hex_tests, HEX)
	elif("chex" in game):
		chex_tests = getTests('chex_tests')
		runTests(chex_tests, CHEX)
	elif("gale" in game):
		gale_tests = getTests('gale_tests')
		runTests(gale_tests, GALE)
	
	

	




	

	

	

		

    

