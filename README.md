# SoftwareCarpentry_LazorProject
This is a lazor project for Software Carpentry 



#code (following is how we design and edit the code)

from itertools import permutations, combinations
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import time


class Lazor():

    '''
    *** CODE PROTOCOL ***
    This code manages all game-related data within "info_dict," a dictionary that stores necessary information for solving the puzzle.

    Functions should be executed in the specified order to ensure proper initialization:
        Lazor = Lazor()
        info_dict = Lazor.read_file('dark_10.bff')
        info_dict = Lazor.load_map(info_dict)
        info_dict = Lazor.find_path(info_dict)
        info_dict = Lazor.solution(info_dict)
        Lazor.save_img(info_dict) or Lazor.save_txt(info_dict)

    Note: Except for read_file(), all functions require only "info_dict" as input, simplifying debugging and execution.

    *** FUNCTION OVERVIEW ***
        ** read_file **
            Reads a .bff file and generates a dictionary with initial map data, block details, laser configurations, and target points.

        ** load_map **
            Processes initial data from info_dict, creating additional structures like map points, fixed block positions, map dimensions, and initializing empty lists and dictionaries for further use.

        ** location **
            A helper function for solution() and redundant(). Determines block placement at a given point and updates the laser’s direction if a block is reflective.

        ** find_path **
            Calculates the path for each laser based on block positions, storing them as lists of position tuples in a dictionary. It also generates a list of possible block positions along each path for block placement.

        ** generate_all_possible_situations **
            A helper function for solution() (Method 1).
            Inputs:
                - possible_list: List of initial possible block placements.
                - block_list: Ordered list of blocks to be placed, where different orders affect laser paths.
            Outputs:
                Updated possible_list with new block combinations along the laser path.

        ** solution **
            The main function to solve the puzzle using three methods, each suited to different scenarios:
                - Scenario 1: All blocks are used, and they must align along laser paths.
                - Scenario 2: Some blocks remain after reaching all targets, and they should not block laser paths.
                - Scenario 3: Cases like the "dark" series, where blocking certain laser paths is essential.

            Solution methods:
                Method 1: Generates combinations of block positions using generate_all_possible_situations(), ensuring all blocks align along paths.
                Method 2: Removes points beyond the last target on each laser path and identifies potential blank positions. If enough blank spots are available, these placements are used.
                Method 3: Brute-force combination generation to find a valid solution.

        ** possible_comb **
            Helper function for solution(). Initializes info_dict and checks each block placement combination to see if all target points are hit, returning a list of unmet targets if any.

        ** redundant **
            A helper function for solution(), handling Method 2. Returns remaining block positions as a dictionary or False if not enough positions are available.

        ** remove_duplicated_element **
            Utility function to remove duplicate items from a list, used throughout the code.

        ** save_txt **
        ** set_color **
        ** pixel_values_list **
        ** pixel_values_dictionary **
        ** save_img **
    '''

    def read_file(self, filename):
        '''
        Reads a .bff file and generates a dictionary with four key components:

        *** map ***
        A list of lists representing the layout of the map.

        *** block ***
        A dictionary indicating the number of each type of block available.

        *** lazor ***
        A dictionary with laser information, including starting positions and directions.
        Format: {(starting position): [direction]}

        *** target_point ***
        A list of tuples indicating the positions of target points on the map.
        '''

    def load_map(self, info_dict):
        '''
        The load_map() function processes the initial data from read_file(),
        which provides the map layout in a structure like [['o'], ['A']].

        This function also initializes various empty lists and dictionaries
        required for subsequent functions.

        The output is a dictionary containing key information:

        *** map_points ***
            A list of tuples representing positions that lasers can pass through.

        *** fixed_block_position ***
            A dictionary with fixed block locations (e.g., "A", "B", "C"),
            in the format: {(x, y): "type"}.

        *** map_length ***
            The width of the map, stored as an integer.

        *** map_height ***
            The height of the map, stored as an integer.
        '''

        

    def location(
        self, position_x, position_y, direction_x, direction_y
    ):
        '''
        This function is separated from find_path() to simplify control flow
        and reduce nested if-else conditions.

        It determines that at the top or bottom edge of a block, the x-direction
        changes, while at the left or right edge, the y-direction changes.
    '''


    def find_path(self, info_dict):
        '''
                The find_path() function generates the path for each laser and identifies
                potential positions for placing additional blocks.

                *** Output ***
                    ** info_dict['find_path'] **
                        A dictionary where each laser path is stored as a list of position tuples,
                        labeled by the starting point of each laser.

                    ** info_dict['lazor'] **
                        A dictionary storing the starting points and directions of the original lasers
                        as well as any generated by refractive blocks.

                    ** info_dict['passed_blocks'] **
                        A dictionary tracking the blocks each laser passes through, used to identify
                        viable block positions, as some blocks may be hit multiple times.

                    ** info_dict['possible_block_position'] **
                        A list of locations along the laser path where additional blocks can be placed.

                *** Method ***
                    ** Path Calculation **
                        - Starting from each laser’s initial position, uses reflect_block_position() to check
                          if a laser hits a block at each step.
                            - If a block is hit, the laser changes direction or stops (if it hits a type B block).
                            - If no block is hit, the point is added to the possible block list.
                        - Adds each point along the path to info_dict['find_path'].
                        - Continues in the direction of the laser until it goes out of range.

                    ** Refractive Blocks **
                        - Lasers created by refractive blocks (type C) are stored in a separate dictionary,
                          c_lazor, and added to info_dict['lazor'] for rerunning the function with the new lasers.
                        - Ensures that duplicate lasers are not added to c_lazor.

                    ** Starting Point Check **
                        - Handles cases where lasers are blocked on both sides by "A" or "B" blocks at their starting point
                          to prevent confusion in path calculations.

                    ** Infinite Loop Prevention **
                        - Prevents infinite loops caused by laser paths trapped between four blocks by setting a counter,
                          which breaks the loop if it exceeds 100 iterations.
        '''


    def generate_all_possible_situations(self, info_dict, possible_list, block_list):
        '''
                The generate_all_possible_situations() function generates updated lists of possible
                block placements for the next block in sequence.
        '''

    def solution(self, info_dict):
        '''
                *** SOLUTION METHODS ***

                    ** METHOD 1 **
                        - Situation 1: All blocks are used, and each block must align with a laser path.

                        - Generate permutations of blocks based on order.
                        - Start with the initial map and explore possible placements for the first block.
                        - For each combination of block placements, check if all target points are hit.
                        - If successful:
                            - Check for remaining blocks. If any remain, proceed to Method 2; otherwise, return info_dict.
                        - If unsuccessful:
                            - Use generate_all_possible_situations() to place the next block, then repeat.
                        - If no combination of block types satisfies the conditions, proceed to Method 3.

                    ** METHOD 2 **
                        - Situation 2: Some blocks remain after all target points are hit; remaining blocks must not block laser paths.

                        - Refer to redundant() to handle this scenario.

                    ** METHOD 3 **
                        - Situation 3: Specific scenarios like "dark" series puzzles, where some laser paths must be blocked.
                          This method also serves as a fallback if previous methods do not solve the problem.

                        - Use itertools.combinations to generate all possible placements on blank positions.
                        - Use remove_duplicated_elements() to eliminate combinations with multiple blocks in the same position.
                        - Loop through and evaluate all combinations to find a solution.
        '''


    def possible_comb(self, info_dict, p, comb):
        '''
        This function generates laser paths based on the given list of possible positions (p)
        and the specified block order (comb).

        Returns a list indicating the number of target points not intersected by any laser path.
        '''


    def redundant(self, info_dict):
        '''
        This function updates info_dict by generating a 'redundant_list' of
        positions where remaining blocks cannot be placed. It follows a similar
        approach as in the solution() function when generating possible positions.

        *** METHOD ***
            *   Move all existing block locations.
            *   Check for any remaining blocks, addressing scenarios that combine
                both method 1 and method 2 from solution().
            *   Remove points beyond the last target point on each laser path,
                allowing block placement only before these points.
            *   Use all_possible_situation() to create a list of positions that
                are unsuitable for block placement.
            *   Exclude unsuitable positions from the blank positions list. If the
                available positions are fewer than the remaining blocks, return False;
                otherwise, return a dictionary of valid positions for remaining blocks.
        '''

    def remove_duplicated_element(self, listA):
        '''
        Small utility function to remove duplicate elements from a list.
        '''
        
    def save_txt(self, info_dict):
        '''
        Generates a .txt file based on the data in info_dict,
        which is produced by the solution function.
        '''

    def change_to_pixel_value(self, a, blockSize=50, gapSize=5):
        '''
        Converts a given value to its corresponding pixel value.
        This function is frequently used in subsequent functions.
        '''


    def pixel_values_list(self, info_dict):
        '''
        This function is to generate a list of pixel values of target points.
        '''

    def pixel_values_dictionary(self, info_dict):
        '''
        Generates a dictionary of pixel values for the laser paths.

        *** info_dict ***
        A dictionary created by previous functions, containing information
        such as filename, map layout, lasers, and target points.
        This data is used to create an image showing the solution to the 'lazor' game.
        '''


    def save_img(self, info_dict):
        '''
        Generates a .png image file to visually display the solution to the game.

        *** info_dict ***
        A dictionary generated by previous functions, containing details such as
        filename, map layout, lasers, and target points. This data is used to
        create an image that illustrates how the 'lazor' game is solved.
        '''




if __name__ == "__main__":
    start_time = time.time()  # Start timing

    a = Lazor()
    b = a.read_file('yarn_5.bff')
    b = a.load_map(b)
    b = a.find_path(b)
    c = (a.solution(b))
    a.save_txt(c)
    a.save_img(c)

    end_time = time.time()  # End timing
    print(f"Total execution time: {end_time - start_time:.2f} seconds")






