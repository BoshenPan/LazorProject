from itertools import permutations, combinations
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import time


class Lazor():

    def read_file(self, filename):
        '''
        This function reads a .bff file, which contains the setup information for
        a Lazor game level, and stores it in a dictionary format. The dictionary
        will contain the map layout, block details, lazor configurations, and
        target points needed for the game.

        Parameters:
            filename (str): The path to the .bff file to read.

        Returns:
            file (dict): A dictionary containing the parsed game setup information.
        '''
        # Open the specified file.
        f = open(filename, "r")
        file = {}  # Initialize an empty dictionary to store game information.
        file['filename'] = filename[:-4]  # Remove file extension and store the name.

        # Initialize map as a list of lists that will store the board layout.
        map = []

        # Skip lines until "GRID START" is found.
        for line in f:
            if line.strip("\n") == "GRID START":
                break  # Exit loop once "GRID START" line is reached.

        # Read the map layout between "GRID START" and "GRID STOP".
        for line in f:
            if line.strip("\n") == "GRID STOP":
                break  # Exit loop once "GRID STOP" line is reached.

            # Create a list for each line in the map.
            lst = []
            for x in line.strip("\n"):
                if x != " ":
                    lst.append(x)
            map.append(lst)  # Add each row of the map to the main map list.

        # Add the map layout to the dictionary.
        file["map"] = map

        # Initialize dictionaries for block and lazor information and a list for target points.
        block = {}
        lazor = {}
        point = []

        # Process remaining lines to capture block, lazor, and target point data.
        for line in f:
            # Check for lines defining blocks of types A, B, or C.
            if line.startswith('A') or line.startswith('B') or line.startswith('C'):
                lst = line.strip("\n").split(" ")
                block[lst[0]] = int(lst[1])

            # Check for lines defining lazors, marked by 'L'.
            if line.startswith('L'):
                lst = line.strip("\n").split(" ")
                lazor[(int(lst[1]), int(lst[2]))] = [int(lst[3]), int(lst[4])]
                # Lazor dictionary stores start position (x, y) as key and direction [vx, vy] as value.

            # Check for lines defining target points, marked by 'P'.
            if line.startswith('P'):
                lst = line.strip("\n").split(" ")
                point.append((int(lst[1]), int(lst[2])))  # Append each target point as a tuple (x, y).

        # Add block, lazor, and target point information to the dictionary.
        file["block"] = block
        file["original_lazor"] = lazor
        file["target_point"] = point

        # Close the file after reading all content.
        f.close()

        return file  # Return the dictionary with all game setup information.

    def load_map(self, info_dict):
        '''
        This function processes the initial map information and generates
        additional details required for solving the Lazor puzzle. It updates
        the info_dict dictionary with the map dimensions, possible laser path
        points, fixed block positions, and blank positions where blocks can be placed.

        Parameters:
            info_dict (dict): A dictionary containing initial map information
                              loaded from a .bff file.

        Returns:
            info_dict (dict): Updated dictionary with additional map information.
        '''

        # Initialize lists and dictionaries in info_dict to store various positions and paths.
        info_dict['map_points'] = []  # All points where the laser can pass.
        info_dict['blank_position'] = []  # Positions where blocks can be placed.
        info_dict['fixed_block_position'] = {}  # Positions of fixed blocks (types A, B, C).
        info_dict['find_path'] = {}  # Used later to store paths found for each laser.
        info_dict['block_position'] = {}  # Stores current positions of blocks.
        info_dict['lazor'] = {}  # Lazor start positions and directions.
        info_dict['possible_block_position'] = []  # Valid positions for new blocks.

        # Load the raw map information from the 'map' key.
        raw_info = info_dict["map"]

        # Calculate the map's length and height in grid units (each grid space is split into two for positioning).
        info_dict["map_length"] = 2 * len(raw_info[0])
        info_dict["map_height"] = 2 * len(raw_info)

        # Iterate over each row (y) and column (x) in the map.
        for y in range(len(raw_info)):
            for x in range(len(raw_info[y])):
                # Check if the position does not allow blocks ('x' represents no-block zones).
                if raw_info[y][x] == 'x':
                    pass  # Skip this cell if no blocks are allowed.
                else:
                    # Calculate the four corner points surrounding each grid cell.
                    possible_points = [
                        (2 * x + 1, 2 * y),
                        (2 * x, 2 * y + 1),
                        (2 * x + 2, 2 * y + 1),
                        (2 * x + 1, 2 * y + 2)
                    ]

                    # Add each point to 'map_points' if it’s not already included.
                    for p in possible_points:
                        if p not in info_dict['map_points']:
                            info_dict['map_points'].append(p)

                    # If the cell allows blocks ('o'), mark it as a blank position for possible block placement.
                    if raw_info[y][x] == 'o':
                        info_dict['blank_position'].append(
                            (2 * x + 1, 2 * y + 1)
                        )
                    # If a fixed block (A, B, C) is in the cell, add its position to 'fixed_block_position'.
                    else:
                        info_dict["fixed_block_position"][
                            (2 * x + 1, 2 * y + 1)] = raw_info[y][x]

        # Copy fixed block positions to 'block_position' for current block setup.
        for f in info_dict['fixed_block_position']:
            info_dict['block_position'][f] = \
                info_dict['fixed_block_position'][f]

        # Add target points to 'map_points' to ensure they are recognized as accessible positions.
        for target_point in info_dict['target_point']:
            if target_point not in info_dict['map_points']:
                info_dict['map_points'].append(target_point)

        # Add laser starting points to 'map_points' to ensure they are recognized as accessible positions.
        for lazor_start in info_dict['original_lazor']:
            if lazor_start not in info_dict['map_points']:
                info_dict['map_points'].append(lazor_start)

        # Return the updated dictionary containing all map information.
        return info_dict

    def location(
            self, position_x, position_y, direction_x, direction_y
    ):
        '''
        This function calculates the new location and direction of a laser
        when it hits a reflective block at a given position. The reflection
        behavior depends on whether the laser hits the block horizontally or
        vertically.

        Parameters:
            position_x (int): The x-coordinate of the current laser position.
            position_y (int): The y-coordinate of the current laser position.
            direction_x (int): The x-component of the laser's direction.
            direction_y (int): The y-component of the laser's direction.

        Returns:
            dict: A dictionary where the key is the new position (x, y) of the
                  laser after reflection, and the value is a list representing
                  the new direction [direction_x, direction_y] after reflection.
        '''

        # If the laser hits the block on the top or bottom (horizontal reflection),
        # reverse the x-direction of the laser's movement.
        if position_y % 2 == 1:
            return {(position_x + direction_x, position_y):
                        [direction_x * -1, direction_y]}

        # If the laser hits the block on the left or right (vertical reflection),
        # reverse the y-direction of the laser's movement.
        elif position_x % 2 == 1:
            return {(position_x, position_y + direction_y):
                        [direction_x, direction_y * -1]}

    def find_path(self, info_dict):
        '''
        This function calculates and tracks the paths of lasers on the board,
        updating directions when lasers hit blocks (reflect, opaque, or refract).
        It also handles generating new lasers when a refract block is hit.

        Parameters:
            info_dict (dict): A dictionary containing map information, block
                              positions, lasers, and target points.

        Returns:
            info_dict (dict): Updated dictionary with laser paths, passed blocks,
                              and potential block positions.
        '''

        # Initialize necessary fields in info_dict to store paths and laser interactions.
        info_dict['find_path'] = {}  # Stores the path of each laser.
        c_lazor = {}  # Temporary storage for lasers created by refract blocks.
        info_dict['passed_blocks'] = {}  # Tracks blocks passed by each laser.
        info_dict['lazor'].update(info_dict["original_lazor"])  # Copy original lasers to 'lazor' dictionary.

        # Iterate over each laser's starting position and direction.
        for i in info_dict['lazor']:
            if i not in info_dict['find_path']:
                info_dict['find_path'][i] = []

            # Set initial laser position and direction.
            x, y = i
            direction_x, direction_y = info_dict["lazor"][i]
            info_dict['passed_blocks'][i] = []

            times = 0  # Counter to prevent infinite loops if lasers get trapped.

            # Continue tracing the laser path until it moves out of map bounds.
            while 0 <= x <= info_dict['map_length'] and \
                    0 <= y <= info_dict['map_height']:
                # Calculate the next location and direction if the laser hits a reflective block.
                (key, value), = Lazor.location(
                    self, x, y, direction_x, direction_y
                ).items()

                # Check if the next position contains a block.
                if key in info_dict["block_position"]:
                    # Record the block if it's the first time the laser has encountered it.
                    if key not in info_dict['passed_blocks'][i] and \
                            info_dict["block_position"][key] != "B":
                        info_dict['passed_blocks'][i].append(key)
                        info_dict['possible_block_position'] = []

                    # Reflective block (type A): Change direction accordingly.
                    if info_dict["block_position"][key] == "A":
                        direction_x, direction_y = value

                    # Opaque block (type B): Block the laser completely.
                    elif info_dict["block_position"][key] == "B":
                        if (x, y) in info_dict["map_points"]:
                            info_dict['find_path'][i].append((x, y))
                        break  # Stop the loop as the laser cannot pass an opaque block.

                    # Refractive block (type C): Generate a new laser while continuing with the current one.
                    elif info_dict["block_position"][key] == "C":
                        if (x + direction_x, y + direction_y) not in \
                                info_dict["lazor"]:
                            c_lazor[(x + direction_x, y + direction_y)] = \
                                (direction_x, direction_y)
                        # Continue tracing the current laser.
                        direction_x, direction_y = value

                # Add the current position to the laser's path.
                if (x, y) in info_dict["map_points"]:
                    info_dict['find_path'][i].append((x, y))

                # Check if the next position is a viable block position.
                if key not in info_dict['possible_block_position'] and \
                        key in info_dict['blank_position'] and \
                        key not in info_dict["block_position"]:
                    info_dict['possible_block_position'].append(key)

                # Move the laser to the next position.
                x += direction_x
                y += direction_y
                times += 1  # Increment counter to prevent infinite loops.

                # Break the loop if too many steps are taken.
                if times > 100:
                    break

        # Recheck each original laser's starting point to ensure it's not blocked on both sides.
        for rsb in info_dict['original_lazor']:
            position_x, position_y = rsb

            # Check for blocks on the vertical sides of the starting point.
            if (position_x, position_y + 1) in info_dict['block_position'] and \
                    (position_x, position_y - 1) in info_dict['block_position']:
                up = (position_x, position_y - 1)
                down = (position_x, position_y + 1)
                if (info_dict['block_position'][up] == "A" or "B") and \
                        (info_dict['block_position'][down] == "A" or "B"):
                    info_dict['find_path'][rsb] = [rsb]

            # Check for blocks on the horizontal sides of the starting point.
            elif (position_x - 1, position_y) in info_dict['block_position'] and \
                    (position_x + 1, position_y) in info_dict['block_position']:
                left = (position_x - 1, position_y)
                right = (position_x + 1, position_y)
                if (info_dict['block_position'][left] == "A" or "B") and \
                        (info_dict['block_position'][right] == "A" or "B"):
                    info_dict['find_path'][rsb] = [rsb]

        # If no refracted lasers were created, return the updated info_dict.
        if not c_lazor:
            return info_dict
        else:
            # Add new lasers generated by refract blocks to 'lazor' and run the function again.
            for c in c_lazor:
                info_dict['lazor'][c] = c_lazor[c]
            return Lazor.find_path(self, info_dict)

    def generate_all_possible_situations(self, info_dict, possible_list, block_list):
        '''
        This function generates all possible combinations of block placements
        along laser paths, ensuring that each configuration is valid. It updates
        the block positions on the map and checks if the lasers pass through
        each required point.

        Parameters:
            info_dict (dict): Dictionary containing map and game information.
            possible_list (list): List of possible positions for placing blocks.
            block_list (list): List of blocks to be placed.

        Returns:
            new_list (list): A list of new valid block position combinations.
        '''

        # Initialize a list to store valid new combinations of block placements.
        new_list = []

        # Iterate over each set of possible positions in the input list.
        for i in possible_list:
            info_dict['find_path'] = {}
            info_dict['lazor'] = {}
            info_dict['block_position'] = {}

            # Update block positions with fixed blocks and the current block configuration.
            info_dict['block_position'].update(info_dict['fixed_block_position'])
            info_dict['block_position'].update(zip(i, block_list))

            # Generate the path for lasers with the current block configuration.
            new_path = Lazor.find_path(self, info_dict)

            # After placing new blocks, some previously passed blocks may no longer be passed.
            # Collect all passed blocks from the current path configuration.
            passed_blocks = []
            for pb in info_dict['passed_blocks']:
                passed_blocks += info_dict['passed_blocks'][pb]

            # Check if all blocks in the current combination `i` are passed by the lasers.
            judge = [False for block in i if block not in passed_blocks]

            # If there are still possible block positions and all blocks in `i` are passed:
            if new_path['possible_block_position'] != [] and judge == []:
                # Add new positions to the list for further exploration.
                for p in new_path['possible_block_position']:
                    new_list.append(i + [p])  # Append the new block position to the combination.
            else:
                # If `possible_block_position` is empty, the combination is invalid
                # because there's no space for the next block.
                pass

        # Return the list of valid new combinations for further exploration.
        return new_list

    def solution(self, info_dict):
        '''
        This function attempts to find a solution for the Lazor puzzle by trying
        various block placements on the board. It uses multiple methods to check
        different configurations of moveable blocks until all target points are hit.

        Parameters:
            info_dict (dict): Dictionary containing the board setup, block positions,
                              laser paths, and other necessary information.

        Returns:
            info_dict (dict): Updated dictionary containing the final block placements
                              if a solution is found.
        '''

        # Initialize the block positions and list of moveable blocks.
        info_dict['block_position'] = {}
        moveable_blocks = []

        for i in info_dict['block']:
            if i != 'B':
                moveable_blocks += [i] * info_dict['block'][i]

        # Generate all unique permutations of the moveable blocks.
        arranged_moveable_blocks = list(
            set(permutations(moveable_blocks, len(moveable_blocks))))

        # Check each combination of block placements.
        while arranged_moveable_blocks != []:
            comb = arranged_moveable_blocks[0]
            i = 0

            info_dict['find_path'] = {}
            info_dict['lazor'] = {}
            info_dict['block_position'] = {}
            info_dict['block_position'].update(info_dict['fixed_block_position'])
            info_dict = Lazor.find_path(self, info_dict)

            # Store initial possible block positions as lists for generate_all_possible_situations().
            possible_list = [[x] for x in info_dict['possible_block_position']]

            # Try placing each block in the combination.
            while i + 1 < len(comb):
                possible_list = Lazor.generate_all_possible_situations(
                    self, info_dict, possible_list, comb)

                # Check each generated possible position configuration.
                for p in possible_list:
                    # Check if the current combination hits all target points.
                    judge = Lazor.possible_comb(self, info_dict, p, comb)

                    # If all target points are hit, check for redundant blocks.
                    if judge == []:
                        rb = Lazor.redundant(self, info_dict)

                        # If there are no redundant blocks (Method 1 successful).
                        if rb is False:
                            pass

                        # If redundant blocks exist, apply Method 2.
                        else:
                            info_dict['block_position'].update(rb)
                            info_dict = Lazor.find_path(self, info_dict)

                            path = []
                            for p in info_dict['find_path']:
                                path += info_dict['find_path'][p]

                            # Verify that all target points are in the path.
                            judge = [
                                False for c in info_dict['target_point']
                                if c not in path
                            ]
                            # Return info_dict if all targets are hit.
                            if judge == []:
                                return info_dict
                            else:
                                pass
                i += 1
            # Remove the current combination from the list and proceed with the next.
            arranged_moveable_blocks.remove(comb)

        # Method 3: When other methods fail.
        moveable_blocks = []
        able_positions = []
        violent_combinations = []

        # Gather all positions where blocks can be placed.
        for bp in info_dict['blank_position']:
            able_positions.append(bp)

        # Generate combinations of blocks in all able positions.
        for block in info_dict['block']:
            moveable_blocks += [block] * info_dict['block'][block]

            # Generate initial combinations if violent_combinations is empty.
            if violent_combinations == []:
                violent_combinations += list(
                    set(combinations(able_positions, info_dict['block'][block])))
            else:
                # Add new combinations to existing ones if there are previous combinations.
                new_list = []
                for vc in violent_combinations:
                    append = list(
                        combinations(able_positions, info_dict['block'][block]))
                    for ap in append:
                        new_list.append(vc + ap)
                violent_combinations = new_list

        listed_violent_combinations = []
        for c in violent_combinations:
            test = list(c)
            test = Lazor.remove_duplicated_element(self, test)
            if len(test) == len(moveable_blocks):
                listed_violent_combinations.append(test)

        # Check each remaining combination to find a valid solution.
        for remain in listed_violent_combinations:
            judge = Lazor.possible_comb(self, info_dict, remain, moveable_blocks)
            if judge == []:
                return info_dict

    def possible_comb(self, info_dict, p, comb):
        '''
        This function checks whether a given combination of block positions
        (p) and block types (comb) successfully hits all target points.
        It updates the paths based on the current block configuration and
        verifies if all targets are covered.

        Parameters:
            info_dict (dict): Dictionary containing the board setup, block
                              positions, lasers, and other necessary information.
            p (list): List of positions where blocks are placed.
            comb (list): List of block types placed at each position in `p`.

        Returns:
            judge (list): An empty list if all target points are covered;
                          otherwise, a list with `False` values for any missed
                          target points.
        '''

        # Initialize paths and laser configurations in info_dict for the current combination.
        info_dict['find_path'] = {}
        info_dict['lazor'] = {}
        info_dict['block_position'] = {}

        # Update block positions with the current combination of positions and block types.
        info_dict['block_position'].update(zip(p, comb))

        # Add fixed block positions to the block position dictionary.
        info_dict['block_position'].update(info_dict['fixed_block_position'])

        # Calculate the paths of lasers based on the current block configuration.
        info_dict = Lazor.find_path(self, info_dict)

        # Collect all points in the laser paths into a single list.
        path = []
        for p in info_dict['find_path']:
            path += info_dict['find_path'][p]

        # Remove any duplicate points in the path to avoid redundant checks.
        path = Lazor.remove_duplicated_element(self, path)

        # Check if all target points are covered in the path.
        # If any target point is not covered, add `False` to `judge`.
        judge = [
            False for c in info_dict['target_point']
            if c not in path
        ]

        # Return `judge`. If empty, all targets are hit; otherwise, some are missed.
        return judge

    def redundant(self, info_dict):
        '''
        This function checks for redundant block placements and identifies
        positions that cannot be used for placing remaining blocks. It analyzes
        laser paths and target points to determine whether enough valid positions
        are available for all moveable blocks.

        Parameters:
            info_dict (dict): Dictionary containing the board setup, blank
                              positions, block positions, and laser paths.

        Returns:
            dict or bool: If there are enough valid positions, returns a dictionary
                          with positions for remaining blocks. If not, returns False.
        '''

        # Initialize lists for moveable blocks and positions where blocks can be placed.
        moveable_blocks = []
        able_positions = []

        # Add all blank positions where blocks can be placed.
        for bp in info_dict['blank_position']:
            able_positions.append(bp)

        # Populate the list of moveable blocks based on available block quantities.
        for block in info_dict['block']:
            moveable_blocks += [block] * info_dict['block'][block]

        # Remove fixed block positions from the list of moveable blocks and able positions.
        for fixed_block in info_dict['block_position']:
            if fixed_block not in info_dict['fixed_block_position']:
                moveable_blocks.remove(info_dict['block_position'][fixed_block])
                able_positions.remove(fixed_block)

        # Method 1: If no moveable blocks remain, return an empty dictionary.
        if moveable_blocks == []:
            return {}

        # Method 2: If there are remaining blocks, check if there are enough valid positions.
        else:
            unable_positions = []

            # Analyze each laser path to find target points and remove unnecessary positions.
            for lazor in info_dict['find_path']:
                test = info_dict['find_path'][lazor].pop()

                # Remove points after the last target point in the path.
                while test not in info_dict['target_point'] and \
                        info_dict['find_path'][lazor] != []:
                    test = info_dict['find_path'][lazor].pop()

                if info_dict['find_path'][lazor] == []:
                    pass
                else:
                    # Recalculate the path from the laser starting point.
                    x, y = lazor[0], lazor[1]
                    direction_x, direction_y = info_dict['lazor'][lazor]

                    while (x, y) in info_dict['find_path'][lazor]:
                        (key, value), = Lazor.location(
                            self, x, y, direction_x, direction_y
                        ).items()

                        if key not in info_dict['block_position'] and \
                                key in info_dict['blank_position']:
                            unable_positions.append(key)
                        else:
                            # Update the direction if a block was encountered.
                            direction_x, direction_y = value
                        x += direction_x
                        y += direction_y

            # Remove duplicate entries in unable_positions.
            unable_positions = Lazor.remove_duplicated_element(self, unable_positions)

            for up in unable_positions:
                able_positions.remove(up)

            # Check if there are enough able positions for the remaining blocks.
            if len(able_positions) < len(moveable_blocks):
                return False  # Not enough positions to place all blocks.
            else:
                # Return a dictionary mapping able positions to remaining moveable blocks.
                return zip(able_positions, moveable_blocks)

    def remove_duplicated_element(self, listA):

        return sorted(set(listA), key=listA.index)

    def save_txt(self, info_dict):
        '''
        This function saves the current state of the Lazor board, including
        the positions of blocks, to a text file. The output file provides a
        visual representation of the board layout with block placements.

        Parameters:
            info_dict (dict): Dictionary containing board information, including
                              the map layout and block positions.

        Returns:
            file (file object): The file object for the saved text file.
        '''

        # Set the filename using the base name from info_dict with '.txt' extension.
        filename = info_dict['filename'] + '.txt'

        f = open(filename, 'w')

        map = info_dict['map']
        block_position = info_dict['block_position']

        # Iterate over each row (y) in the map.
        for y in range(len(map)):
            # Iterate over each column (x) in the row.
            for x in range(len(map[0])):
                # Check if the current position has a block placed.
                # Blocks are placed at positions (2*x+1, 2*y+1) to match grid coordinates.
                if (2 * x + 1, 2 * y + 1) in block_position:
                    # Write the block type (e.g., 'A', 'B', 'C') to the file.
                    f.write(block_position[(2 * x + 1, 2 * y + 1)] + ' ')
                else:
                    # Write the original map cell value ('o', 'x') if no block is placed.
                    f.write(map[y][x] + ' ')

            # Write a newline after each row to format the board correctly.
            f.write('\n')

        f.close()

        return f

    def change_to_pixel_value(self, a, blockSize=50, gapSize=5):
        '''
        This function converts a grid coordinate to a pixel value, taking into
        account the size of each block and the gap between blocks.

        Parameters:
            a (int): The grid coordinate to be converted to pixel value.
            blockSize (int): The size of each block in pixels (default is 50 pixels).
            gapSize (int): The size of the gap between blocks in pixels (default is 5 pixels).

        Returns:
            int: The pixel value corresponding to the grid coordinate.
        '''

        # Calculate the pixel position by accounting for the block and gap sizes.
        a = gapSize + (blockSize + gapSize) * a

        return a

    def pixel_values_list(self, info_dict):

        target_point = info_dict['target_point']
        tp = []
        for (x, y) in target_point:
            x = Lazor.change_to_pixel_value(self, x) / 2
            y = Lazor.change_to_pixel_value(self, y) / 2
            tp.append((x, y))
        return tp

    def pixel_values_dictionary(self, info_dict):

        find_path = info_dict['find_path']
        lp = {}
        for (x0, y0) in find_path:
            x1 = Lazor.change_to_pixel_value(self, x0) / 2
            y1 = Lazor.change_to_pixel_value(self, y0) / 2
            lst = []
            for (x, y) in find_path[(x0, y0)]:
                x = Lazor.change_to_pixel_value(self, x) / 2
                y = Lazor.change_to_pixel_value(self, y) / 2
                lst.append((x, y))
            lp[(x1, y1)] = lst
        return lp

    def save_img(self, info_dict):

        # Generate a new list of the pixel values of map.
        map = info_dict['map']
        block_position = info_dict['block_position']
        for y in range(len(map)):
            for x in range(len(map[0])):
                if (2 * x + 1, 2 * y + 1) in block_position:
                    map[y][x] = block_position[(2 * x + 1, 2 * y + 1)]

        # Create a new image.
        w_blocks = len(map[0])
        h_blocks = len(map)
        SIZE = (Lazor.change_to_pixel_value(self, w_blocks),
                Lazor.change_to_pixel_value(self, h_blocks))
        # Set the background of the image.
        img = Image.new("RGB", SIZE, (99, 90, 91))

        # Show three different types of blocks:
        # A - reflect blocks, B - opaque blocks, C - refract blocks
        # and fixed blocks of these three types.
        for y, row in enumerate(map):
            for x, block_ID in enumerate(row):
                if block_ID != 'x':
                    block = Image.open(block_ID + '.jpg')
                    block = block.resize((50, 50))
                    pxl = block.load()
                    for a in range(50):
                        for b in range(50):
                            img.putpixel((Lazor.change_to_pixel_value(self, x) + a,
                                          Lazor.change_to_pixel_value(self, y) + b),
                                         pxl[a, b])

        # Draw lines to show the path of the laser.
        draw = ImageDraw.Draw(img)
        lazor = info_dict['lazor']
        original_lazor = info_dict['original_lazor']
        find_path = Lazor.pixel_values_dictionary(self, info_dict)

        # For laser touching block A/B.
        for (x, y) in lazor:
            if (x, y) in original_lazor:
                x = Lazor.change_to_pixel_value(self, x) / 2
                y = Lazor.change_to_pixel_value(self, y) / 2
                draw.ellipse((x - 3, y - 3, x + 3, y + 3),
                             fill=(255, 0, 0))
                lst = find_path[(x, y)]
                draw.line(lst, width=5, fill=(255, 0, 0))

        # For laser passing block C.
        for (x, y) in lazor:
            if (x, y) not in info_dict['original_lazor']:
                x = Lazor.change_to_pixel_value(self, x) / 2
                y = Lazor.change_to_pixel_value(self, y) / 2
                lst = find_path[(x, y)]
                direction_x = lst[1][0] - lst[0][0]
                direction_y = lst[1][1] - lst[0][1]
                lst.insert(0, (x - direction_x, y - direction_y))
                draw.line(lst, width=5, fill=(255, 0, 0))

        # Draw points to show the position of target points.
        target_point = Lazor.pixel_values_list(self, info_dict)
        for (x, y) in target_point:
            draw.ellipse((x - 5, y - 5, x + 5, y + 5), fill=(255, 0, 0))
            draw.ellipse((x - 2, y - 2, x + 2, y + 2), fill=(255, 255, 255))
        filename = info_dict['filename']
        img.save(filename + '.png')

        # Display the image with matplotlib
        plt.imshow(img)
        plt.axis('off')  # Turn off axes for a clean image display
        plt.show()


if __name__ == "__main__":
    start_time = time.time()  # Start timing

    a = Lazor()
    b = a.read_file('mad_1.bff')
    b = a.load_map(b)
    b = a.find_path(b)
    c = (a.solution(b))
    a.save_txt(c)
    a.save_img(c)

    end_time = time.time()  # End timing
    print(f"Total execution time: {end_time - start_time:.2f} seconds")






