-------------------------------------------------------------------------------
---------------------------------- GUIDELINES ---------------------------------
-------------------------------------------------------------------------------

local cell_width, cell_height = 10, 10

-------------------------------------------------------------------------------------
-------------------------------- GRID FUNCTIONS -------------------------------------
-------------------------------------------------------------------------------------
-- Contains (1) the cell coordinates and (2) the initial winding number increment
local event_list = {}

local function fixLineWindingNumber(line, start, end)
    -- Loop through line changing winding number 'till the net winding number is zero
    -- RETURN VOID
end

local function computeGridDimension(scene)
    -- this should return a "optimal" width and height after
    return cell_width, cell_height
end

local function makeGrid(window_width, window_height, cell_width, cell_height)
    -- returns an empty grid (a bidimensional table), where each cell 
    -- contains (1) its bounding box and (2) a table with the intersecting segments
    -- and the (3) initial winding number

    -- RETURN: A TABLE WITH FORMAT CELL[i][j] = {xmin, ymin, xmax, ymax, initialWindingNumber, segments = {} }
end

local function intersectSegmentCell(x0, y0, x1, y1, segment)
    -- 1) Check path bounding box against Cell
    --      -> if it is fully inside, then return true
    --      -> if it is not, check if one of the extreme control points are inside the cell.
    -- 2) Ray cast -> check for intersection with right side
    --      -> If it does not intersect, check for intersection with top side (by rotating the cell)

    -- RETURN : BOOLEAN
end

local function walkInPath(path, grid)
    -- For each segment inside the path, walk through it in a Brenseham/Tripod-fashion
    -- push segment into intersecting cell
    -- write events to event_list
    -- (1) pf final control point inside? 
    -- (2) Is the segment going up or down
    -- (3) Test respective cells. If up: test cell[i][j+1], cell[i][j-1], cell[i+1][j]; 
    --      otherwise, cell[i-1][j], cell[i][j+1], cell[i][j-1]
    --          -> if cell[i+1][j], event_list.push(i+1, j, increment)
    --          -> if cell[i-1][j], event_list.push(i-1, j, increment)
    -- (4)  

    -- RETURN: VOID
end

local function prepareGrid(scene)
    -- 1) Compute grid dimensions
    -- 2) Create grid
    -- 3) loop through paths inside scene
    -- 4) Sort event_list -> insertion_sort (or any other stable sort)
    -- 5) fix initial winding numbers

    -- RETURN: FILLED GRID
end

local function getCell(x, y, grid)
    -- Compute coordinates of cell containing (x,y)
    -- RETURN: i, j -> INTEGERS!!!
end

local function export_cell(i,j,cell)
    -- 4) compose scene
    -- 5) Export svg

    -- RETURN: VOID, but exports a .svg file
end