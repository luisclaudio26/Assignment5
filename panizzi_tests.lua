--[[
	Alana Rasador Panizzi
	panizzis@gmail.com
	Rio de Janeiro, 15 jan 2016 
	]]

	math.randomseed(31415)

	local driver = require"driver"
	local image = require"image"
	local chronos = require"chronos"
	local unpack = table.unpack
	local floor = math.floor
	local max = math.max
	local min = math.min
	local abs = math.abs
	local _M = driver.new()
	local xform = require"xform"
	local bezier = require"bezier"
	local util = require"util"
	local significant = util.significant
	local distinct = util.distinct
	local similar = util.similar
	local xform_meta = _M.meta
	local cut2 = bezier.cut2
	local cut3 = bezier.cut3
	local cut2rc = bezier.cut2rc
	local canonic2r = bezier.canonic2r
	local split2r = bezier.split2r
	local bernstein = require"bernstein"
	local quadr = require"quadratic"
	local quadratic = quadr.quadratic
	local sin = math.sin
	local asin = math.asin
	local sqrt = math.sqrt
	local ceil = math.ceil
	local pi = math.pi
	local atan = math.atan
	local cub = require"cubic"
	local cubic = cub.cubic
	local concat = table.concat
	local initialtangent = bezier.initialtangent3
	local finaltangent = bezier.finaltangent3
	local blue = require"blue"

-- OMG this' so awful
local GRID

-- Fifth Assignment

	-------------------------------------------------------------------------------------
	-------------------------------- GRID FUNCTIONS -------------------------------------
	-------------------------------------------------------------------------------------l]

	-- Contains (1) the cell coordinates and (2) the initial winding number increment
	local event_list = {}
	local cell_size = 10

	--DONE AND WORKING
	local function fixLineWindingNumber(cells, event_list)
		local str = ""
		local w = event_list[#event_list][1]
		local y, x = event_list[#event_list][2], event_list[#event_list][3]
		for i = #event_list-1, 1, -1 do
			next_y = event_list[i][2]
			if next_y == y then
				next_x = event_list[i][3]
			else
				next_x = 1
			end
			if w ~= 0 then
				for dx = next_x + 1, x-1 do
					cells[y][dx].initialWindingNumber = w

					-- This is a "fake" segment, just to correctly
					-- paint the cell
					local fake_seg = { foo = function(a, b, c) return 0 end }
					cells[y][dx].shapes = { {["segment"] = fake_seg, ["path"] = event_list[i][4]} }

					str = str.."["..dx.."]["..y.."].w = "..w..", "
				end
			else
				--print("w=0")
			end
			if next_y ~= y then 
				w = event_list[i][1]
				--print(str)
				str = ""
			else 
				w = w + event_list[i][1] 
			end
			y, x = event_list[i][3], next_y
		end
		if w ~= 0 then
			for dx = 1, x do
				cells[y][dx].initialWindingNumber = w

				local fake_seg = { foo = function(a, b, c) return 0 end }
				cells[y][dx].shapes = { {["segment"] = fake_seg, ["path"] = event_list[i][4]} }

				str = str.."["..dx.."]["..y.."].w = "..w..", "
			end
			--print(str)
		end
		-- Loop through line changing winding number 'tilatl the net winding number is zero
	    -- RETURN VOID
	end

	--DONE AND WORKING
	local function computeGridDimension(scene)
	    -- this should return a "optimal" width and height after
	    -- We will need the viewport info about width and height
	    -- added to render: scene.width, scene.height =  width, height

	    -- Now we will make a fixed grid, with a size cell_size + extra_space
	    -- All cells must have the same size, so..
	    local extra_space_width = (scene.width%cell_size)/floor(scene.width/cell_size)
	    local extra_space_height = (scene.height%cell_size)/floor(scene.height/cell_size)
	    return cell_size+ extra_space_width , cell_size+ extra_space_height
	end

	--DONE AND WORKING
	local function makeGrid(window_width, window_height, cell_width, cell_height)
	    -- returns an empty grid (a bidimensional table), where each cell 
	    -- contains (1) its bounding box and (2) a table with the intersecting segments
	    -- and the (3) initial winding number
	    local grid_width, grid_height = window_width/cell_width, window_height/cell_height
	    local cells = {["width"]  = cell_width, ["height"]  = cell_width}
	    for i = 0, grid_width + 1 do
	    	cells[i] = {}
	    	for j = 0, grid_height + 1 do
	    		cells[i][j] = {
	    		["xmax"] = cell_width*(i), ["xmin"] =  cell_width*(i-1),
    			["ymax"] = cell_height*(j), ["ymin"] = cell_height*(j-1),
    			["initialWindingNumber"] = 0, ["border"] = false, 
    			["shapes"] = {}, } -- [LC] These are not segments, but shapes (segment + paint)

    		end
    		-- Define a border, so we can clip it when needed
    		cells[i][0].border = true
    		cells[i][grid_width + 1].border = true
    	end
    	for j = 0, grid_width + 1 do 
    		cells[0][j].border = true
    		cells[grid_height + 1][j].border = true
    	end
	    -- RETURN: A TABLE WITH FORMAT CELL[i][j] = {xmin, ymin, xmax, ymax, initialWindingNumber, segments = {} }
	    return cells
	end

	--DONE AND WORKING	
	local function findCellCoord(x, y, cells)
		-- Why not ceil()?
		return floor(x/cells.width) + 1, floor(y/cells.height) + 1
	end

	--DONE AND WORKING
	----------------------- uses my definition of segment and rayCastingFunction-----------------------
	local function rayCasting(segment, xmin, xmax, ymin, ymax)
		local rayCastingFunction = segment.foo
		-- Tests if it intersect
		if segment.ymax >= ymax and segment.ymin <= ymax then
			if abs(rayCastingFunction(segment, xmin, ymin) + 
				rayCastingFunction(segment, xmax, ymin))  == 1  or
			   abs(rayCastingFunction(segment, xmin, ymax) +   
				rayCastingFunction(segment, xmax, ymax)) == 1 then
				return true
			end
		end
		return false
	end
	---------------------------------------------------------------------------------------------------

	---------------------- uses my definition of segment and rayCastingFunction----------------
	--DONE AND WORKING
	local function create_mirror_segment(segment)
		local y, x = segment.x, segment.y
		local new_segment = {}
		local foo = segment.foo
		new_segment.x, new_segment.y = x, y
		if foo == winding_number_linear then
			local a, b = y[2] - y[1], x[1] - x[2]
			local c =  -(a*x[1] + b*y[1])
			new_segment.a, new_segment.b, new_segment.c = a, b, c
			new_segment.foo = foo
			new_segment.xmax, new_segment.xmin, 
				new_segment.ymax, new_segment.ymin  = bounding_values(x, y)
			return new_segment
		end
		if foo == winding_implicit_quadratic then
			w = segment.w
			if w then typ = rquadratic else typ = quadratic end
		else -- cubic
			typ = "cubic"
		end
		p = {}
		p.k, p.segments = 0, {}
		monotonic_segment_implicit(p, 0, foo, x, y, typ, w)
		return p.segments[0]
	end
	-------------------------------------------------------------------------------------------------
	--DONE AND WORKING FOR TESTED EXAMPLES (MAYBE WE CAN TRY TO FIND ERRORS LATER)
	local function intersectSegmentCell(i, j, cells, segment)
	    -- 1) Check path bounding box against Cell
	    --      -> if it is fully inside, then return true
	    --      -> if it is not, check if one of the extreme control points are inside the cell.
	    ----------------------- uses my definition of segment -----------------------
	    local x, y = segment.x, segment.y
	    -----------------------------------------------------------------------------
	    local i1, j1 = findCellCoord(x[1], y[1], cells)
	    local ifinal, jfinal = findCellCoord(x[#x], y[#y], cells)
	    -- case final or first control point is inside cell
	    if (i1 == i and j1 == j) or (ifinal == i and jfinal == j) then
	    	return true
	    end

	    -- 2) Ray cast 
		--> check for intersection with right side
	    --> If it does not intersect, check for intersection with top side (by rotating the cell)
	    local xmin, xmax, ymin, ymax = cells[i][j].xmin, cells[i][j].xmax, 
	    	cells[i][j].ymin, cells[i][j].ymax
		if rayCasting(segment, xmin, xmax, ymin, ymax) then
			return true
		end
		inv_segment = create_mirror_segment(segment)  -- Inverse ray casting
		if rayCasting(inv_segment, ymin, ymax, xmin, xmax) then
			return true
		end
		return false
		--RETURN : BOOLEAN
	end

	--DONE AND WORKING
	local function inside_cell(x, y, i, j, cells)
		xi, yi = findCellCoord(x, y, cell)
		return xi == i and yj == j
	end

	----------- uses my definition of segment (what a segment must have) --------------
	-- DONE AND WORKING
	local function createClosingCellSegment(x, y)

		local a, b = y[2] - y[1], x[1] - x[2]
		local c =  -(a*x[1] + b*y[1])

		return {["x"] = x, ["y"] = y, ["foo"] = winding_number_linear, 
			["sign"] = -1, ["a"] = a, ["b"] = b, ["c"] = c}

	end
	-----------------------------------------------------------------------------------

	--DONE AND WORKING
	local function walkInPath(path, cells)
	    -- For each segment inside the path, walk through it in a Brenseham/Tripod-fashion
	    -- push segment into intersecting cell
	    -- write events to event_list

	    for l, segment in ipairs(path.shape.segments) do
	    	----------------------- uses my definition of segment -----------------------
	    	local x, y = segment.x, segment.y

	    	-----------------------------------------------------------------------------
	    	local finali, finalj = findCellCoord(x[#x], y[#y], cells)
	    	-- First visited cell
	    	local i, j = findCellCoord(x[1], y[1], cells)
	    	-- treat case we start outside viewport
	    	message1 = "other"
	    	--(2) Is the segment going up or down, lef or right
	    	if (finalj >= j) then segment_going_up = true
	    		else segment_going_up = false end
	    	if (finali >= i) then segment_going_right = true
	    		else segment_going_right = false end
	    	if segment_going_up then
		    	if segment_going_right then
		    		message = ("This segment is going up right")
		    	else
		    		message = ("This segment is going up left")
		    	end
		    elseif segment_going_right then
	    		message = ("This segment is going down right")
	    	else
    			message = ("This segment is going down left")
    		end

		   	while true do
		    	table.insert(cells[i][j].shapes, {["segment"] = segment, ["path"] = path} )
		    	print("Inserting in cell ", i, j, segment.type, segment.x[1], segment.y[1], segment.x[2], segment.y[2])
		    	-- cell that contains the first control point
	    	-- (1) pf final control point inside this cell or did we reach a border of the viewport?
		    	if (i == finali and j == finalj) or cells[i][j].border then 
		    		break --leaves the while loop and goes to another segment
		    	end
		    -- (3) Test respective cells. If up: test cell[i][j+1], cell[i][j-1], cell[i+1][j]; 
	    		if segment_going_up then 
	    			if intersectSegmentCell(i, j+1, cells, segment) then
	    				flag = true
			    	--	print("caseup")
		    			-- leave above, increment initial winding number
		    			table.insert(event_list, {1, i, j, path})
	    				--print(i, j, "->", i, j+1)
		    			i, j = i, j+1
		    		end
		    	elseif intersectSegmentCell(i, j-1, cells, segment) then -- leave under, decrement initial winding number  
	    			flag = true
		    	--	print("casedown")
	    			table.insert(event_list, {-1, i, j, path})
	    			--print(i, j, "->", i, j-1)
	    			i, j = i, j-1
    			end 
    			--if not flag then
		    		if segment_going_right then
			    		if intersectSegmentCell(i+1, j, cells, segment) then
				    	--	print("caseright")		
		    				flag = true
			    			-- leave through the right, store a new closing segment if segment's end is not over the cell
			    			if  cells[i][j].ymax > y[#y] then

			    				-- [LC] I think this was wrong! x1 = x2 = x[#x], 
			    				-- and y1 = y[#y], y2 = cells[i][j].ymax
			    				local newx, newy = {x[#x], x[#x]}, {y[#y], cells[i][j].ymax}

			    				local s = createClosingCellSegment(newx, newy)
			    				table.insert(cells[i][j].shapes, {["segment"] = s, ["path"] = path})
			    			end
			    			--print(i, j, "->", i+1, j)
		    				i, j = i+1, j
			    		end
			    	elseif intersectSegmentCell(i-1, j, cells, segment) then
			    		--print("caseleft")
		    			--leave through the left, do nothing
		    			flag = true
		    			--print(i, j, "->", i-1, j)
						i, j = i-1, j
			    	end
		    	--end
		    	if not flag then 
		    		print("errrooooooooooo")
		    		break
		    	end
		    	--print("new cell", i, j)
	    	end
	    end
	        -- RETURN: VOID
	end
	
	--DONE AND WORKING
	--needs a key as well so we order this array by this particular table position
	local function CountingSort(array, key)
		local count = {}
	    for j = 1, #array do
	    	i = array[j][key]
	    	if not count[i] then count[i] = 0 end
	        count[i] = count[i] + 1
	    end
	    local total = 1
		for i  = 1, #array do
			if not count[i] then count[i] = 0 end
		    local oldCount = count[i]
		    count[i] = total
		    total = total + oldCount
		end
		local output = {}
		for i, x in pairs(array) do
		    output[count[x[key]]] = x
		    count[x[key]] = count[x[key]] + 1
		end
		return output
	end

	--DONE AND WORKING
	local function prepareGrid(scene)
	    -- 1) Compute grid dimensions
	    local cell_width, cell_height = computeGridDimension(scene)
	    -- 2) Create grid
		cells = makeGrid(scene.width, scene.height, cell_width, cell_height)
	    -- 3) loop through paths inside scene
		-- insert here yout respective sample loop
		for i, e in ipairs(scene.inverse_shapes_list) do
			if e.shape.type == "path" then
				print("enter path")

				-- [LC] Tem de ter as informações do paint também!
				walkInPath(e, cells)
			else
				-- we have to think about it
				-- maybe a function walkInCircle/Triangle/Polygon or segment them 
			end
		end
	    -- 4) Sort event_list -> insertion_sort (or any other stable sort)
		local x_sorted_event_list = CountingSort(event_list, 2) -- (sort by x line)
		local y_sorted_event_list = CountingSort(event_list, 3) -- (sort by x line)
		--[[for i in pairs(y_sorted_event_list) do
			print(i, y_sorted_event_list[i][1], y_sorted_event_list[i][2], y_sorted_event_list[i][3] )
		end--]]

		-- 5) fix initial winding numbers
		fixLineWindingNumber(cells, y_sorted_event_list)	
	    -- RETURN: FILLED GRID

	    return cells
	end


	--TODO
	local function export_cell(i,j,cell)
	    -- 4) compose scene
	    -- 5) Export svg

	    -- RETURN: VOID, but exports a .svg file
	end

-- First Assignment

	function applies(xv, yv, xf, wv)
		local x, y, w = {}, {}, nil
		for i, v in ipairs(xv) do
			x[i], y[i] = xf:apply(xv[i], yv[i])
			if i == 2 and wv then
				x[i], y[i], w = xf:apply(xv[i], yv[i], wv)
			end
		end
		return x, y, w
	end

	function bounding_values(x, y)
		local xmax, xmin, ymax, ymin = x[1], x[1], y[1], y[1]
		for i, v in ipairs(x) do
			xmax, xmin = max(x[i], xmax), min(x[i], xmin)
			ymax, ymin = max(y[i], ymax), min(y[i], ymin)
		end 
		return xmax, xmin, ymax, ymin 
	end

	function colinear(x, y)
		if (not x[4]) then
			return (x[1]*(y[2]- y[3]) + x[2]*(y[3] - y[1]) + x[3]*(y[1] - y[2]) == 0)
		else
			return
			(x[1]*(y[2]- y[3]) + x[2]*(y[3] - y[1]) + x[3]*(y[1] - y[2]) == 0) and
			(x[4]*(y[2]- y[3]) + x[2]*(y[3] - y[4]) + x[3]*(y[4] - y[2]) == 0)
		end
	end

	local function prepare_circle(scene, c)
		local tr, sc, center
	    tr, sc  = xform.translate(c.cx, c.cy), xform.scale(c.r)
		center = tr*sc
		x, y, xbox, ybox = {}, {}, {1, -1, -1, 1}, {1, -1, 1, -1}
		x, y = applies(xbox, ybox, scene.xf*c.xf*center)
		c.xmax, c.xmin, c.ymax, c.ymin = bounding_values(x, y)
	    assert(scene.xf:det() ~= 0, "invalid transformation")
	    assert(c.xf:det() ~= 0, "invalid transformation")
		c.inv_xf = center:inverse()*c.xf:inverse()*scene.xf:inverse()
		c.inside_shape = inside_circle
	end

	function inside_circle(c, px, py)
		local x, y = c.inv_xf:apply(px, py)
		return (x*x + y*y < 1)
	end

	local function prepare_triangle(scene, t)
		local x, y = {t.x1, t.x2, t.x3}, {t.y1, t.y2, t.y3}
		x, y = applies(x, y, scene.xf*t.xf)
		if(colinear(x, y)) then
		 	print "Warning! Colinear dots"
		end
		t.x, t.y = x, y
		t.xmax, t.xmin, t.ymax, t.ymin = bounding_values(x, y)
		t.inside_shape = inside_triangle
	end

	function inside_triangle(t, px, py)
		-- using the baricentric test
		local d, na, nb, s
		local x, y = t.x, t.y
		x21, x31, qx = (x[2] - x[1]), (x[3] - x[1]), (px - x[1])
		y21, y31, qy = (y[2] - y[1]), (y[3] - y[1]), (py - y[1])
		d  = det(x21, x31, y21, y31)
		na = det(qx, x31, qy, y31)
		nb = det(x21, qx, y21, qy)
		s = d/math.abs(d)
		return (s*na > 0 and s*nb > 0 and s*(na + nb) <= s*d)
	end

	local function prepare_polygon(scene, p)
		local i, k, a, b
		local x, y, newpath, coord, path = {}, {}, {}, {}, p.data
		local xf = p.xf*scene.xf
		for i in pairs(path) do
			if (i%2 == 0) then
				k = i/2
				x[k], y[k] = xf:apply(path[i-1], path[i])
				if (i >= 4) then
					a, b = y[k] - y[k-1], x[k-1] - x[k]
					c =  -(a*x[k-1] + b*y[k-1]) 
					newpath[i/2 -1] = {["x"] = {x[k-1], x[k]}, ["y"] = {y[k-1], y[k]}, 
					["a"] = a, ["b"] = b, ["c"] = c }				
				end
			end
		end
		a, b = y[1] - y[k], x[k] - x[1]
		c =  -(a*x[k] + b*y[k]) 
		newpath[0] = {["x"] = {x[k], x[1]}, ["y"] = {y[k], y[1]}, 
			["a"] = a, ["b"] = b, ["c"] = c	}
		p.path = newpath
		p.xmax, p.xmin, p.ymax, p.ymin = bounding_values(x, y)
		p.inside_shape = inside_polygon
		p.a, p.b = y[2] - y[1], x[1] - x[2]
		p.c = -(a*x[1] + b*y[1]) 
	end 

	function inside_polygon(p, px, py, fill_type)
		local result, sign = 0, 0
		for i, e in pairs(p.path) do
			sign = winding_number_linear(e, px, py)
			result = result + sign
		end
		if fill_type == "eofill" then
			return (math.abs(result) % 2 ~= 0)
		elseif fill_type == "fill" then 
			return result ~= 0
		end
		return true
	end

	function winding_number_linear(e, px, py)
		local c, s, bool, expression
		local x, y, sign = e.x, e.y, 0 
		local a, b, c = e.a, e.b, e.c

		--print(x[1], y[1], x[2], y[2], a, b, c)

		if (a < 0 or (a == 0 and b > 0)) then
			bool = y[2] < py and py <= y[1]
		else
			bool = y[1] < py and py <= y[2]
		end
		if bool then	
			if     a ~= 0 then s = util.sign(a)
			elseif b ~= 0 then s = util.sign(b)
			else return 0 end
			a = s*a 	b = s*b 	c = s*c
			if a*px + b*py + c <= 0 then
				sign = s
			end
		end	
		return sign
	end

-- Second Assignment


	-- calculates bounding values in a path
	local function path_bounding_values(p, xmax, ymax, xmin, ymin)
		if( not p.xmax) then 
			return xmax, ymax, xmin, ymin
		else
			return max(xmax, p.xmax), max(ymax, p.ymax), min(xmin, p.xmin), 
				min(ymin, p.ymin)
		end
	end

	function prepare_degenerate(i, p, x, y)
	end

	function prepare_LS(i, j, scene, p, vx, vy, bool)

		local x, y, foo, sign, xmax, xmin, ymax, ymin

		-- [LC] could be replace for vy = vy or blabla
		if vy == nil then
			vx, vy = {p.data[j], p.data[j+2]}, {p.data[j+1], p.data[j+3]}
		end
		if not bool then
			x, y = applies(vx, vy, scene.xf*p.xf)
		else
			x, y = vx, vy
		end

		xmax, xmin, ymax, ymin = bounding_values(x, y)
		foo = winding_number_linear
		if (y[2] - y[1] >= 0) then
			sign = 1
		elseif (y[2] - y[1] < 0 ) then
			sign = -1
		end
		local a, b = y[2] - y[1], x[1] - x[2]
		local c =  -(a*x[1] + b*y[1])
		p.segments[i + p.k] = {["x"] = x, ["y"] = y,
			["xmin"] = xmin, ["xmax"] = xmax, ["ymax"] = ymax, ["ymin"] = ymin, 
			["foo"] = foo, ["sign"] = sign,
			["a"] = a, ["b"] = b, ["c"] = c, ["sign"] = sign, ["type"] = "linear_segment" }

		-- [LC] Add TYPE info to segment - but it seems a good idea to do a blackbox implementation
		-- in which is not necessary to know the type

		p.xmax, p.ymax, p.xmin, p.ymin = 
		path_bounding_values(p, xmax, ymax, xmin, ymin)
	end 

	local function create_segment(p, i, foo, u, v, r)
		umax, umin, vmax, vmin = bounding_values(u, v)
		local sign = 0
		p.xmax , p.ymax, p.xmin, p.ymin = 
			path_bounding_values(p, umax, vmax, umin, vmin)
		if (v[#v] - v[1] >= 0) then
			sign =  1
		elseif (v[#v] - v[1] < 0 ) then
			sign = -1
		end

		p.segments[i+p.k] = {["x"] = u, ["y"] = v, 
			["xmax"] = umax, ["xmin"] = umin,
			["ymax"] = vmax, ["ymin"] = vmin, 
			["foo"] = foo, ["sign"] = sign, ["i"] = i+p.k,
		}
		p.segments[i+p.k].rootx = {}
		p.segments[i+p.k].b = u[1] - u[#u]
		if r then
			p.segments[i+p.k].w = r
		end	
		return p.segments[i+p.k]
	end

	function prepare_QS(i, j, scene, p)
		local num_roots, x, y
		local roots, a, b = {}, {}, {}
		local foo = winding_number_quadratic
		x, y = {p.data[j]  , p.data[j+2], p.data[j+4]}, 
			   {p.data[j+1], p.data[j+3], p.data[j+5]} 
		if(colinear(x, y)) then
			prepare_LS(i, j, scene, p, {x[1], x[3]}, {y[1], y[3]})
			return
		end
		x, y = applies(x, y, scene.xf*p.xf)
		a = {x[1] - 2*x[2] + x[3], y[1] - 2*y[2] + y[3]}
		b = {-(x[2] - x[1]),  -(y[2] - y[1])}
		num_roots = 0
		for l = 1, 2 do
			if(a[l] ~= 0 and b[l]/a[l] >= 0 and b[l]/a[l] <= 1) then
				num_roots = num_roots + 1
				roots[num_roots] = b[l]/a[l]
			end
		end
		table.sort(roots)
		prev_root = 0
		if (num_roots == 0) then
			create_segment(p, i, foo, x, y)
		else 
			roots[0], roots[num_roots + 1] = 0, 1
			for l = 1, num_roots do
				if (roots[l] < 1 and roots[l] > 0 and prev_root ~= roots[l]) then
					u0, v0, u1, v1, u2, v2 = 
						cut2(prev_root, roots[l], x[1], y[1], x[2], y[2], x[3], y[3])
					create_segment(p, i, foo, {u0, u1, u2}, {v0, v1, v2})
					prev_root = roots[l]
					p.k = p.k + 1
				end
			end
			u0, v0, u1, v1, u2, v2 = 
				cut2(prev_root, 1, x[1], y[1], x[2], y[2], x[3], y[3])
			s = create_segment(p, i, foo, {u0, u1, u2}, {v0, v1, v2})
			p.segments[i+p.k].ay, p.segments[i+p.k].by, p.segments[i+p.k].cy = 
				v0 - 2*v1 + v2, -2*v0 + 2*v1, v0
			p.segments[i+p.k].ax, p.segments[i+p.k].bx, p.segments[i+p.k].cx = 
				u0 - 2*u1 + u2, -2*u0 + 2*u1, u0 
			prev_root = roots[l]
		end
	end

	function degenerated_test(e, x, y)

		return e.y == y and e.x == x
	end

	function prepare_CS(i, j, scene, p)
		foo = winding_number_cubic
		x, y = {p.data[j]  , p.data[j+2], p.data[j+4], p.data[j+6]}, 
			   {p.data[j+1], p.data[j+3], p.data[j+5], p.data[j+7]} 
		if (colinear(x, y)) then
			prepare_LS(i, j, scene, p, {x[1], x[4]}, {y[1], y[4]})
			return
		end
		x, y = applies(x, y, scene.xf*p.xf)
		xmin, ymin, xmax, ymax = bounding_values(x, y)
		ax, ay =  3*(x[2] - x[1]) - 6*(x[3] - x[2]) + 3*(x[4] - x[3]), 
		          3*(y[2] - y[1]) - 6*(y[3] - y[2]) + 3*(y[4] - y[3])
		bx, by = -6*(x[2] - x[1]) + 6*(x[3] - x[2]), 
		         -6*(y[2] - y[1]) + 6*(y[3] - y[2])
		cx, cy =  3*(x[2] - x[1]), 
			      3*(y[2] - y[1])
		t = {} s = {}
		roots = {}
		rootx, t[1], s[1], t[2], s[2] = quadratic(ax, bx, cx)	
		rooty, t[3], s[3], t[4], s[4] = quadratic(ay, by, cy)
		k = 0
		for l = 1, 4 do
			if (t[l]) then
				if (s[l] ~= 0) then
					k = k+1
					roots[k] = t[l]/s[l]
				end
			end
		end
		table.sort(roots)
		prev_root = 0
		nx, ny =  x, y
		if (k == 0) then
			create_segment(p, i, foo, x, y)
		else 
			for l = 1, k do
				if (roots[l] < 1 and roots[l] > 0 and roots[l] ~= prev_root) then
					u1, v1, u2, v2, u3, v3, u4, v4  = cut3(prev_root, roots[l], 
						x[1], y[1], x[2], y[2], x[3], y[3], x[4], y[4])
					create_segment(p, i, foo, {u1, u2, u3, u4}, {v1, v2, v3, v4})
					prev_root = roots[l]
						p.k = p.k + 1	
				end
			end
			u1, v1, u2, v2, u3, v3, u4, v4  = cut3(prev_root, 1, 
				x[1], y[1], x[2], y[2], x[3], y[3], x[4], y[4])
			create_segment(p, i, foo, {u1, u2, u3, u4}, {v1, v2, v3, v4})
		end	
	end

	function prepare_RQS(i, j, scene, p)
		foo = winding_number_rational
		x0, y0 = scene.xf:apply(p.xf:apply(p.data[j], p.data[j+1]))
		x1, y1, w = scene.xf:apply(p.xf:apply(p.data[j+2], p.data[j+3], p.data[j+4]))
		x2, y2 = scene.xf:apply(p.xf:apply(p.data[j+5], p.data[j+6]))
		assert(w > 0 and p.data[j+4], "W must be greater than zero")
		if (colinear({x0, x1, x2}, {y0, y1, y2})) then
			prepare_LS(i, j, scene, p, {x0, x2}, {y0, y2}, true)
			return
		end
			roots = {}
			ax, ay = x0 - 2*x1 + x2, y0 - 2*y1 + y2
			bx, by =  -(x1 - x0),  -(y1 - y0)
			k = 0
			if(ax ~= 0 and bx/ax > 0 and bx/ax < 1) then
				k = k+1
				roots[k] = bx/ax
			end
			if(ay ~= 0 and by/ay > 0 and by/ay < 1) then
				k = k+1	
				roots[k] = by/ay
			end
			table.sort(roots)
			if (k == 0) then
				create_segment(p, i, foo, {x0, x1, x2}, {y0, y1, y2}, w)
			else 
				roots[0] = 0
				roots[k+1] = 1
				for l = 1, k do
					if (roots[l] < 1 and roots[l] > 0) then
						u1, v1, r1, u2, v2, r2, u3, v3, r3 = 
						split2r(roots[l], x0, y0, 1, x1, y1, w, x2, y2, 1)
						x0, y0, u0, v0, r, u1, v1 = 
						canonic2r(x0, y0, 1, u1, v1, r1, u2, v2, r2)
						create_segment(p, i, foo, {x0, u0, u1}, {y0, v0, v1},  r )
						p.k = p.k + 1
						u2, v2, u3, v3, r, x2, y2 = 
						canonic2r(u2, v2, r2, u3, v3, r3, x2, y2, 1)
						create_segment(p, i, foo, {u2, u3, x2}, {v2, v3, y2}, r)
					end
				end
			end
	end


	function winding_number_quadratic(e, x, y)
		--remembers the respective x for value y of line
		if (not e.rootx[y]) then
			local x0, y0, x1, y1, x2, y2 = 
			e.x[1], e.y[1], e.x[2], e.y[2], e.x[3], e.y[3]
			local ay, by, cy = y0 - 2*y1 + y2, -2*y0 + 2*y1, y0
			xnum_root, t1, s1, t2, s2 = quadratic(ay, by, (cy - y))
			if (num_root == 0) then
				e.rootx[y] = 0
				return 0
			end
			local r1, r2 = t1/s1 , t2/s2
			--idea for corners: count twice the intersection (<=, >=)
			if (r1 >= 0 and r1 <= 1) then
				r = r1
			elseif (r2 >= 0 and r2  <= 1) then
				r = r2
			else
				e.rootx[y] = 0
				return 0
			end
			local ax, bx, cx = x0 - 2*x1 + x2, -2*x0 + 2*x1, x0 
			xr = ax*r*r + bx*r + cx
			e.rootx[y] = xr
		else 
			xr = e.rootx[y]
		end
		if (x <= xr) then
			return e.sign
		end
		e.rootx[y] = 0
		return 0
	end

	function winding_number_rational(e, x, y)
		--remembers the respective x for value y of line
		if (not e.rootx[y]) then
			if x <= e.xmin then return e.sign end
			local x0, y0, x1, y1, w, x2, y2 = e.x[1], e.y[1], e.x[2], e.y[2], e.w, e.x[3], e.y[3]
			local as, bs = y2 - y0, x0 - x2
			local a, b, c = y0 - 2*y1 + y2, -2*y0 + 2*y1, y0
			local wa, wb, wc =  1 - 2*w + 1,  -2 + 2*w, 1
			a, b, c = a - wa*y, b - wb*y, c - wc*y 
			local num_root, t1, s1, t2, s2 = quadratic(a, b, c)
			if (num_root ~= 0) then
				local r1, r2 = t1/s1 , t2/s2 
				if (r1 > 0 and r1 < 1) then
					r = r1
				elseif (r2 > 0 and r2  < 1) then
					r = r2
				else
					e.rootx[y] = 0
					return 0 
				end
			else
				e.rootx[y] = 0
				return 0
			end 
			a, b, c = x0 - 2*x1 + x2, -2*x0 + 2*x1, x0
			wa, wb, wc =  1 - 2*w + 1,  -2 + 2*w,    1
			
			local wr = wa*r*r + wb*r + wc
			xr = (a*r*r + b*r + c)/wr
			e.rootx[y] = xr
		else
			xr = e.rootx[y]
		end
		if x < xr then
			return e.sign
		end
		e.rootx[y] = 0
		return 0
	end

	function winding_number_cubic(e, x, y)
		--remembers the respective x for value y of line
		if (not e.rootx[y]) then
			local x0, y0, x1, y1, x2, y2, x3, y3 =
				e.x[1], e.y[1], e.x[2], e.y[2], e.x[3], e.y[3], e.x[4], e.y[4]
			if x < e.xmin then return e.sign end
			a =   -y0 + 3*y1 - 3*y2 + y3
			b =  3*y0 - 6*y1 + 3*y2
			c = -3*y0 + 3*y1
			d =    y0 - y
			t = {} rt = {}
			s = {}	
			num_root, t[1], s[1], t[2], s[2], t[3], s[3] = cubic(a, b, c, d)
			for l = 1, num_root do
				rt[l] = t[l]/s[l]
				if (rt[l] >= 0 and rt[l] <= 1) then
					r = rt[l]
					break
				end
			end
			if (not r) then
				e.rootx[y] = 0
				return 0 
			end
			a =   -x0 + 3*x1 - 3*x2 + x3
			b =  3*x0 - 6*x1 + 3*x2
			c = -3*x0 + 3*x1
			d =    x0
			xr = a*r^3 + b*r^2 + c*r + d
			e.rootx[y] = xr
		else
			xr =  e.rootx[y]
		end
		if (x <= xr) then
			return e.sign
		end
		e.rootx[y] = 0
		return 0
	end

-- Third Assignment

	function extract_ramp_colors(ramp, opacity)
		local paints, t, last_color = {}, {}, 0
		for i, e in ipairs(ramp) do
			if i%2 == 1 then  
				t[(i+1)/2] =  e
			else	
				paints[i/2] = e
				paints[i/2][4] = paints[i/2][4]*opacity
				last_color = i/2
			end
		end
		return paints, t, last_color
	end

	function create_array_ramp(size, paints, t, last_color)
		local ramp, k = {}, 1
		--dt is the distance between each color in the ramp
		local dt = t[2] - t[1] 
		if t[1] == t[2] then
			dt = 1
		end
		for i =  1, size do
			ti = (i- 1)/size
			ramp[i] = adds(mult(paints[k],   1 - ((ti - t[k])/dt)), 
			    mult(paints[k+1], 1 - ((t[k+1] - ti)/dt)) )
			--print(ramp[i][1], ramp[i][2], ramp[i][3], ramp[i][4])
			last = k+1
			if (ti >= t[k+1] and last_color ~= k+1) then
				i = i - 1
				k = k+1
				dt = t[k+1] - t[k] 
				if t[k] == t[k+1] then
					dt = 1
				end
			end
		end
		ramp[1], ramp[size] = paints[1], paints[last_color]
		return ramp
	end

	function prepare_gradient_linear(scene, paint)
		local data = paint.data
		-- apply xf to dots in the ramp
		local x, y = applies({data.p1[1], data.p2[1]}, 
			{data.p1[2], data.p2[2]}, paint.xf*scene.xf) 
		local dy, dx = y[2] - y[1], x[2] - x[1]
		-- create xf that transforms the ramp into a segment [0,1]
		local h = sqrt(dx^2+dy^2)
		local rot, sc = xform.rotate(dx/h, -dy/h), xform.scale(1/(h))
		local transl = xform.translate(-x[1], -y[1])
		local paints, t, last_color = extract_ramp_colors(data.ramp, paint.opacity)
		local size = 256
		paint.data.xfnew = sc*rot*transl
		paint.data.c = create_array_ramp(size, paints, t, last_color)
		paint.foo, paint.data.t, paint.data.size = coloring_linear, t, size
	end

	function prepare_gradient_radial(scene, paint)
		local cx, cy, fx, fy, dx, dy, rot, sc, paints, t, k, transl, ramp
		local data, xf = paint.data, paint.xf
		local xfscene = scene.xf
		cx, cy = data.center[1], data.center[2]
		fx, fy = data.focus[1], data.focus[2]
		radius = data.radius
		-- create xf that transforms in unitary circle, with focus in [0,1]
		dx, dy = fx - cx, fy - cy
		inv_scene, sc = scene.xf:inverse() , xform.scale(1/radius)
		transl, inv_xf = xform.translate(-cx, -cy), paint.xf:inverse()
		h = sqrt(dx^2+dy^2)
		if (h ~= 0) then 
			rot = xform.rotate(dx/h, -dy/h)
			paint.data.xfnew = sc*rot*transl*inv_xf*inv_scene
		else
			paint.data.xfnew = sc*transl*inv_xf*inv_scene
		end
		paints, t, last_color = extract_ramp_colors(data.ramp, paint.opacity)
		ramp = {}
		local dq, total_angle, dist_cf = {},  180, sqrt(dx^2 + dy^2)/radius
		for angle = 0, total_angle do 
			if (dist_cf == 0) then dq[angle] = 1
			else dq[angle] = find_q(dist_cf, math.rad(angle))
			end
			size = ceil(dq[angle]*256)
			ramp[angle] = {}
			ramp[angle] = create_array_ramp(size, paints, t, last_color)		
		end 
		paint.data.total_angle, paint.data.dist_cf = total_angle, dist_cf
		paint.data.dq, paint.foo = dq, coloring_radial
		paint.data.t, paint.data.cr, paint.data.size = t, ramp, size
	end

	-- given an angle, finds distance from center to point q that is colinear with
	-- p, dot with an angle alpha with the focus
	function find_q(dx, alpha) 
		local alpha = pi - alpha
		local gamma = asin(dx*sin(alpha))
		if (gamma == 0) then
			return (1 + dx)
		end
		local beta = pi - gamma - alpha
		if (beta == 0) then	return (1 - dx)	end
		return  (sin(beta) / sin(alpha))
	end

	function coloring_solid(data)
		return data
	end

	function spread(x, d, size, spread_type, approx)
		if not approx then approx = ceil end
		if spread_type == "pad" then
			if (x > 0) then 
				return min(approx(x*size), size)
			else 
				return 1
			end
		elseif spread_type == "repeat" then
			return approx(x*d.size % size)
		elseif spread_type == "reflect" then
			if floor(x) % 2 == 0 then
				return approx(x*size % size)
			else
				return approx(d.size - (x*size % size))
			end
		elseif spread_type == "transparent" then
			if 0 < x and x <= 1 then
				return approx(x*size)
			else
				return 1
			end
		end
		print("Spread not supported. Repeat instead")
		return approx(x*d.size % d.size)
	end

	function coloring_linear(data, px, py)
		local d, xf, c = data, data.xfnew, data.c
		local x, y = xf:apply(px, py)
		return c[spread(x, d, d.size, d.ramp.spread)]
	end

	function coloring_radial(data, px, py)
		local d, xf, dq, dist_cf, c, total_angle = data, data.xfnew, data.dq, data.dist_cf, data.cr, data.total_angle
		local x, y = xf:apply(px, py)
		local dx, dy = x - dist_cf, y
		if dx == 0 then angle = 90 else
			angle = ceil(math.deg(atan2(dy,dx)))
			if angle < 0 then print("what") angle = abs(angle) end
		end
		local t = sqrt(dx^2 + dy^2)/dq[angle]
		return c[angle][spread(t, d, ceil(dq[angle]*256), d.ramp.spread)]
	end

	function in_bounding_box(shape, x, y)
		if not (shape.xmin and shape.ymin and shape.xmax and shape.ymax) then
			return true
		end
		return (x >= shape.xmin and x <= shape.xmax and y >= shape.ymin and 
			y <= shape.ymax)
	end

	local function sample(scene, x, y) 
		local colors, k, c =  {}, 1, {}

		-- Get corresponding cell
		local cell_x, cell_y = findCellCoord(x, y, GRID)
		local cell = GRID[cell_x][cell_y]

		for i, e in ipairs(cell.shapes) do
			local seg = e.segment

			local winding_number = seg.foo(seg, x, y) + cell.initialWindingNumber

			if winding_number ~= 0 then
				c = e.path.paint.foo(e.path.paint.data, x, y, e.paint)
				colors[k] = c
				k = k+1
				if (c[4] == 1) then
					break 
				end
			end
		end
		colors[k] = {1, 1, 1, 1}
	    return compose_colors(colors, k)
	end

	function compose_colors(c, k)
		local transp, colors, color = 0, c, {1, 1, 1, 1}
		if (colors[1][4] == 1) then 
			return colors[1]
		end
		if (k >= 2) then
			color = colors[k]
			for i = (k-1), 1, -1  do
				transp = colors[i][4]
				color = adds(mult(colors[i], transp), mult(color, (1 - transp)))
			end
		end
		color[4] = 1
		return color
	end

-- Fourth assignment
-- Texture:

	function prepare_texture(scene, paint, shape)
		local img, spread, xf = paint.data.image, paint.data.spread, paint.xf
		local w, h = img.width , img.height
		paint.data.w, paint.data.h = w, h
		local c = {}
		for i = 1, w do
			c[i] = {}
			for j = 1, h do
				local r, g, b, a = 0, 0, 0, 1
				r, g, b, a = img:get(i, j) 
				c[i][j] = {r, g, b, a}
			end
		end
		newxf = (scene.xf*xf):inverse()
		paint.data.c = c
		local x1, y1 = newxf:apply(0, 1)
		local x2, y2 = newxf:apply(0, 0)
		local improve = false
		if ((x1 - x2)^2 + (y1 - y2)^2) > 16 then
			improve = false
		end
		paint.data.improve = improve
		paint.data.newxf = newxf
		paint.data.size = w
	end

	function coloring_texture(data, px, py, paint)
		local d, c = data, data.c
		local newxf = paint.data.newxf
		local x, y = newxf:apply(px, py)
		local ix2 = spread(x, d, d.w, d.spread, math.ceil)
		local iy2 = spread(y, d, d.w, d.spread, math.ceil)
		if (data.improve) then 
			local ix1 = spread(x, d, d.w, d.spread, math.floor)
			local iy1 = spread(y, d, d.w, d.spread, math.floor)
			color11, color12 = c[ix1][iy1], c[ix1][iy2] 
			color21, color22 = c[ix2][iy1], c[ix2][iy2] 
			color1 = adds(mult(color11, y%1), mult(color12, 1 - y%1))
			color2 = adds(mult(color21, y%1), mult(color22, 1 - y%1))
			return adds(mult(color1, x%1), mult(color2, 1 - x%1))
		end
		return c[ix2][iy2]
	end

-- Implicitization
	local function prepare_path(scene, p)
		p.segments, p.inside_shape, p.k = {}, inside_path, 0
		for i, e in pairs(p.instructions) do
			x, y, j = {}, {}, p.offsets[i]
			if e == "begin_closed_contour" or  e ==  "begin_open_contour" then
				local len = p.data[j] 
				local l =  p.offsets[i+len]
				local x, y = {p.data[l], p.data[j+1]},  
					   {p.data[l + 1], p.data[j+2]}
				if not (x[1] == x[2] and y[1] == y[2]) then
					prepare_LS(i - p.k, j, scene, p, x, y)	
				else 	p.k = p.k - 1
				end
			elseif e == "linear_segment" then
				prepare_LS(i, j, scene, p)
			elseif e == "quadratic_segment" then 
				prepare_implicit_quadratic(i, j, scene, p, true)
			elseif e == "rational_quadratic_segment" then 
				prepare_implicit_quadratic(i, j, scene, p, false)
			elseif e == "cubic_segment" then 
				prepare_implicit_cubic(i, j, scene, p)
			elseif e == "linear_segment_with_length" then
				local x, y = {p.data[j], p.data[j+2]}, 
					{p.data[j+1], p.data[j+3]}
				local len = p.data[j+4]
				local size = sqrt((x[1] - x[2])^2 + (y[1] - y[2])^2)
				if size ~= 0 and len ~= 0 then
					x[2] = (x[2] - x[1])*len/size + x[1]
					y[2] = (y[2] - y[1])*len/size + y[1]
					prepare_LS(i, j, scene, p, x, y)
				end
			end
		end
		--[[for i, e in pairs(p.segments) do 
			print(i)
			--print("Xmax, xmin, ymax, ymin", e.xmax, e.xmin, e.ymax, e.ymin)
			print("x, y", e.x[1], e.y[1], e.x[2], e.y[2], e.x[3], e.y[3], e.x[4], e.y[4])
		end--]]
	end

	function prepare_implicit_quadratic(i, j, scene, p, integral, xv, yv)
		local x, y, pts, klm, b2p
		x, y, px, py = {}, {}, {}, {}
		if(integral) then
			--print("CASE: Integral" )
			x, y = {p.data[j]  , p.data[j+2], p.data[j+4]}, 
				   {p.data[j+1], p.data[j+3], p.data[j+5]}
			if(colinear(x, y)) then
				prepare_LS(i, j, scene, p, {x[1], x[3]}, {y[1], y[3]})
				return
			end
			x, y, w = applies(x, y, scene.xf*p.xf)
		else 
			--print("CASE: Racional" )
			px, py, pw = {p.data[j]  , p.data[j+2], p.data[j+5]}, 
				         {p.data[j+1], p.data[j+3], p.data[j+6]},
				          p.data[j+4]
			x, y, w = applies(px, py, scene.xf*p.xf, pw)
			if (colinear({x[1], x[2], x[3]}, {y[1], y[2], y[3]})) then
				prepare_LS(i, j, scene, p, {x[1], x[3]}, {y[1], y[3]}, true)
				return
			end
		end
		if xv then x, y = xv, yv end
		roots, num_roots = find_roots_quadratic(x, y)
		local foo = winding_implicit_quadratic
		local typ = "rquadratic"
		if integral then typ = "quadratic" w = nil end
		separate_monotonic_segments(p, i, x, y, roots, num_roots, foo, typ, w)
	end 

	function find_roots_quadratic(x, y)
		local roots, a, b, num_roots = {}, {}, {}, 0
		a = {x[1] - 2*x[2] + x[3], y[1] - 2*y[2] + y[3]}
		b = {-(x[2] - x[1]),  -(y[2] - y[1])}
		for l = 1, 2 do
			if(a[l] ~= 0 and b[l]/a[l] >= 0 and b[l]/a[l] <= 1) then
				num_roots = num_roots + 1
				roots[num_roots] = b[l]/a[l]
			end
		end
		roots[0], roots[num_roots + 1] = 0, 1
		table.sort(roots)
		return roots, num_roots
	end

	function separate_monotonic_segments(p, i, x, y, rt, num_roots, foo, typ, w, klm)
		local t = typ
		local prev_root = 0
		local roots = rt
		if (num_roots == 0) then
			monotonic_segment_implicit(p, i, foo, x, y, t, w)
		else 
			for l = 1, num_roots do
				--print("root", roots[l])
				if (roots[l] < 1 and roots[l] > 0 and prev_root ~= roots[l]) then
					u, v, r = cut(prev_root, roots[l], x, y, t, w)
					monotonic_segment_implicit(p, i, foo, u, v, t, r)
					prev_root = roots[l]
					p.k = p.k + 1
				end
			end
			u, v, r = cut(prev_root, 1, x, y, t, w)
			monotonic_segment_implicit(p, i, foo, u, v, t, r)
		end
	end

	function cut(ta, tb, x, y, typ, w)
		if typ == "quadratic" then
			u0, v0, u1, v1, u2, v2 = 
				cut2(ta, tb, x[1], y[1], x[2], y[2], x[3], y[3])
			return {u0, u1, u2}, {v0, v1, v2}
		elseif typ == "rquadratic" then
			u0, v0, u1, v1, r, u2, v2 = 
				cut2rc(ta, tb, x[1], y[1], x[2], y[2], w, x[3], y[3])
			return {u0, u1, u2}, {v0, v1, v2}, r
		elseif typ == "cubic" then
			u1, v1, u2, v2, u3, v3, u4, v4  = 
			cut3(ta, tb, x[1], y[1], x[2], y[2], x[3], y[3], x[4], y[4])
			return  {u1, u2, u3, u4}, {v1, v2, v3, v4}
		end
	end

	function find_triangle_r(x, y, right)
		xmax, ymax = max(x[1], x[3]) , max(y[1], y[3])
		xmin = min(x[1], x[3])
		if right then
			xt = xmax
			if xmax == x[1] then
				if ymax == y[3] then yt = y[3]
				else yt = y[1] end
			elseif ymax == y[3] then yt = y[3]
				else yt = y[1]
			end
		else
			if xmin == x[1] then
				if ymax == y[3] then yt = y[3]
				else yt = y[1] end
			elseif ymax == y[3] then yt = y[3]
				else yt = y[1]
			end
			xt = xmin
		end
		return {x[1], xt, x[3]}, {y[1], yt, y[3]}
	end

	function adjugate(a, b, c, d, e, f, g, h, i)
		return xform.xform(-f*h + e*i, c*h - b*i, -c*e + b*f, 
				f*g - d*i, -c*g + a*i,  c*d - a*f, 
				-e*g + d*h, b*g - a*h, -b*d + a*e)
	end

	function monotonic_segment_implicit(p, i, foo, u, v, typ, w)
		local index = create_segment(p, i, foo, u, v, w)
		--quadratic precomputations
		if typ == "quadratic" or typ == "rquadratic" then
			if (not w) then w = 1 end
			local pts, klm, b2p = {}, {}, {}
			pts = adjugate(u[1], u[2], u[3], v[1], v[2], v[3], 1, w, 1)
			b2p = xform.xform(0, 0.5, 0, 0, 0, 1, 1, 0, 0)
			local sign = util.sign(v[3] - v[1])
			local a, b = (v[3] - v[1]), (u[1] - u[3])
			local c =  -sign*(a*u[1] + b*v[1])
			a, b = a*sign, b*sign
			p.segments[i+p.k].klm =  b2p*pts
			p.segments[i+p.k].sign = sign
			p.segments[i+p.k].signx = util.sign(u[3] - u[1])
			p.segments[i+p.k].a = {a, b, c}
			if a*u[2] + b*v[2] + c*w > 0 then
				p.segments[i+p.k].triangletoright = true
				--print("dir")
			else
				--print("esq")
				p.segments[i+p.k].triangletoright = false
			end
			if typ == "rquadratic" then 
				local xt, yt = find_triangle_r(u,v, p.segments[i+p.k].triangletoright)
				p.segments[i+p.k].triangle = {["x"] = xt, ["y"] = yt}
			else
				p.segments[i+p.k].triangle = {["x"] = u, ["y"] = v}
			end
		--cubic precomputations
		elseif typ == "cubic" then
			local t = find_bounding_triangle(u, v)
			local sign = util.sign(v[4] - v[1])
			local a, b = (v[4] - v[1]), (u[1] - u[4])
			local c =  -sign*(a*u[1] + b*v[1])
			a, b = a*sign, b*sign
			p.segments[i+p.k].triangle = t
			p.segments[i+p.k].sign = sign
			p.segments[i+p.k].signx = util.sign(u[4] - u[1])
			p.segments[i+p.k].a = {a, b, c}
			if a*t.x[2] + b*t.y[2] + c > 0 then
				p.segments[i+p.k].triangletoright = true
				--print("dir")
			else
				p.segments[i+p.k].triangletoright = false
				--print("esq")
			end
			m = {} m[1] = {} m[2] = {} m[3] = {}
			m[1], m[2], m[3], m[4] = {2, -3, -2}, {-6, 3, 3}, {-2, -3, -2}
			local px1, px2, px3, px4 = u[1], -3*u[1] + 3*u[2], 3*u[1] - 6*u[2] + 3*u[3],
				-u[1] + 3*u[2] - 3*u[3] + u[4]
			local py1, py2, py3, py4 = v[1], -3*v[1] + 3*v[2], 3*v[1] - 6*v[2] + 3*v[3], 
				-v[1] + 3*v[2] - 3*v[3] + v[4]
			px, py  = {px1, px2, px3, px4},
				{py1, py2, py3, py4}
			local midx = (max(u[1], u[4]) + min(u[1], u[4]))/2
			local midy = (max(v[1], v[4]) + min(v[1], v[4]))/2
			coef = prepare_resultant(px, py)
			p.segments[i+p.k].intersectionsign = util.sign(implicit_resultant(coef, midx, midy))
			p.segments[i+p.k].resultant_coefficients = prepare_resultant(px, py)
		end
	end

	function prepare_implicit_cubic(i, j, scene, p, integral)
		local x, y, pts, klm, b2p, roots, num_roots
		x, y = {p.data[j]  , p.data[j+2], p.data[j+4], p.data[j+6]}, 
			   {p.data[j+1], p.data[j+3], p.data[j+5], p.data[j+7]}
		if(colinear(x, y)) then
			prepare_LS(i, j, scene, p, {x[1], x[4]}, {y[1], y[4]})
			return
		end
		x, y = applies(x, y, scene.xf*p.xf)
		if (-x[1] + 3*x[2] - 3*x[3] + x[4] == 0) and (-y[1] + 3*y[2] - 3*y[3] + y[4]==0) then
			xaux = (3/2)*(x[2] - x[1]) + x[1]
			yaux = (3/2)*(y[2] - y[1]) + y[1]
			prepare_implicit_quadratic(i, j, scene, p, true, 
				{x[1], xaux, x[4]}, {y[1], yaux, y[4]})
			return
		end
		roots, num_roots = find_roots_cubic(x, y)
		inflections, num_intersections = find_inflections(x, y, roots, num_roots)
		table.sort(inflections)
		local foo = winding_implicit_cubic
		local typ = "cubic"
		local w = nil
		separate_monotonic_segments(p, i, x, y, inflections, num_intersections, foo, typ, w, klm)
		
	end 

	function find_roots_cubic(x, y)
		local ax, ay =  3*(x[2] - x[1]) - 6*(x[3] - x[2]) + 3*(x[4] - x[3]), 
		          3*(y[2] - y[1]) - 6*(y[3] - y[2]) + 3*(y[4] - y[3])
		local bx, by = -6*(x[2] - x[1]) + 6*(x[3] - x[2]), 
		         -6*(y[2] - y[1]) + 6*(y[3] - y[2])
		local cx, cy =  3*(x[2] - x[1]), 
			      3*(y[2] - y[1])
		t = {} s = {}
		roots = {}
		rootx, t[1], s[1], t[2], s[2] = quadratic(ax, bx, cx)	
		rooty, t[3], s[3], t[4], s[4] = quadratic(ay, by, cy)
		k = 0
		for l = 1, 4 do
			if (t[l]) then
				if (s[l] ~= 0) then
					k = k+1
					roots[k] = t[l]/s[l]
				end
			end
		end
		table.sort(roots)
		return roots, k
	end

	function cross(v1, v2, v3)
		a, b, c, d = v1[1], v1[2], v1[3], v1[4]
		l, e, f, g = v2[1], v2[2], v2[3], v2[4]
		h, i, j, k = v3[1], v3[2], v3[3], v3[4]
		return {d*f*i - c*g*i - d*e*j + b*g*j + c*e*k - b*f*k,
			   -d*f*h + c*g*h + d*l*j - a*g*j - c*l*k + a*f*k,
			    d*e*h - b*g*h - d*l*i + a*g*i + b*l*k - a*e*k,
			   -c*e*h + b*f*h + c*l*i - a*f*i - b*l*j + a*e*j}
	end

	function find_inflections(x, y, roots, num_roots)
		local num_inflections, t, s = num_roots, {}, {}
		local v1, v2, v3 =  
			{x[1], -3*x[1] + 3*x[2], 3*x[1] - 6*x[2] + 3*x[3], -x[1] + 3*x[2] - 3*x[3] + x[4]},
			{y[1], -3*y[1] + 3*y[2], 3*y[1] - 6*y[2] + 3*y[3], -y[1] + 3*y[2] - 3*y[3] + y[4]},
			{1, 0, 0, 0}
		local d = cross(v1, v2, v3)
		num_inf, t[1], s[1], t[2], s[2] = quadratic(1, -(d[3]/d[2]), d[4]/(3*d[2]))
		--double_point
		num_double_inf, t[3], s[3], t[4], s[4] = quadratic(d[2]^2, -d[2]*d[3], d[3]^2 - d[2]*d[4])
		for k = 1, (num_inf + num_double_inf) do
			l = k + (2 - num_inf)
	  		if(s[l] ~= 0 and  t[l]/s[l] >= 0 and t[l]/s[l] <= 1) then
				num_inflections = num_inflections + 1
				table.insert(roots, t[l]/s[l])
			end
		end
		table.sort(roots)
		return roots, num_inflections  
	end

	function det3new(m) 
		local sum = 0
		for i = 1, 3 do
			sum = sum + ((-1)^(i-1))*m[1][i]*det2(restrict_matrix(m, 1, i))
		end
		return sum
	end

	function det3(m)
		return m[1][1]*m[2][2]*m[3][3] + m[1][2]*m[2][3]*m[3][1] + m[1][3]*m[2][1]*m[3][2]
			 - m[1][3]*m[2][2]*m[3][1] - m[1][2]*m[2][1]*m[3][3] - m[1][1]*m[2][3]*m[3][2]
	end

	function det2(m)
		return m[1][1]*m[2][2] - m[1][2]*m[2][1]
	end

	function restrict_matrix(m, a, b)
		local new_m, k, l = {}, 1, 1
		for i, ei in ipairs(m) do
			if (i ~= a) then
				new_m[k] = {} 
				l = 1
				for j, ej in ipairs(ei) do
					if (j ~= b) then
						new_m[k][l] = m[i][j]
						l = l + 1
					end
				end
				k = k + 1 
			end
		end
		return new_m
	end

	function resultant3(px, py, bool)
		local m = {}
		local f0, f1, f2, f3 = px[1], px[2], px[3], px[4]
		local g0, g1, g2, g3 = py[1], py[2], py[3], py[4]
		m[1] = {f1*g0 - f0*g1, f2*g0 - f0*g2,                 f3*g0 - f0*g3}
		m[2] = {f2*g0 - f0*g2, f3*g0 + f2*g1 - f1*g2 - f0*g3, f3*g1 - f1*g3}
		m[3] = {f3*g0 - f0*g3, f3*g1 - f1*g3,                 f3*g2 - f2*g3}
		return det3(m)--]]
		--[[ local a, b, c, d = px[1], px[2], px[3], px[4]
		local e, f, g, h = py[1], py[2], py[3], py[4]
		return -d^3*e^3 + c*d^2*e^2*f - b*d^2*e*f^2 + a*d^2*f^3 - c^2*d*e^2*g + 
	     2*b*d^2*e^2*g + b*c*d*e*f*g - 3*a*d^2*e*f*g - a*c*d*f^2*g - 
	     b^2*d*e*g^2 + 2*a*c*d*e*g^2 + a*b*d*f*g^2 - a^2*d*g^3 + c^3*e^2*h - 
	     3*b*c*d*e^2*h + 3*a*d^2*e^2*h - b*c^2*e*f*h + 2*b^2*d*e*f*h + 
	     a*c*d*e*f*h + a*c^2*f^2*h - 2*a*b*d*f^2*h + b^2*c*e*g*h - 
	     2*a*c^2*e*g*h - a*b*d*e*g*h - a*b*c*f*g*h + 3*a^2*d*f*g*h + 
	     a^2*c*g^2*h - b^3*e*h^2 + 3*a*b*c*e*h^2 - 3*a^2*d*e*h^2 + 
	     a*b^2*f*h^2 - 2*a^2*c*f*h^2 - a^2*b*g*h^2 + a^3*h^3--]]
		 --]]
	
	end


	function prepare_resultant(px, py )
		local d, c, b, a = px[1], px[2], px[3], px[4]
		local h, g, f, e = py[1], py[2], py[3], py[4]
		local d1, d2, d3, d1h1, d1h2, d2h1, h1, h2, h3
		d3 = -e^3
	    d2 = c*e^2*f- b*e*f^2 + a*f^3 +  2*b*e^2*g - 3*a*e*f*g
	    d1 = b*c*e*f*g - a*c*f^2*g - b^2*e*g^2 + 2*a*c*e*g^2 + a*b*f*g^2 - a^2*g^3 - c^2*e^2*g
	    d1h1 = -3*b*c*e^2 + 2*b^2*e*f +  a*c*e*f - 2*a*b*f^2 - a*b*e*g  + 3*a^2*f*g 
	    d1h2 = -3*a^2*e
	    d2h1 =  3*a*e^2
	    h3 = a^3
	    h2 = -b^3*e + 3*a*b*c*e + a*b^2*f - 2*a^2*c*f - a^2*b*g
	    h1 = c^3*e^2  - b*c^2*e*f  + a*c^2*f^2 + b^2*c*e*g - 2*a*c^2*e*g - a*b*c*f*g + a^2*c*g^2 
	    return {d, h, d1, d2, d3, d1h1, d1h2, d2h1, h1, h2, h3}
	end

	function implicit_resultant(coefients, x, y)
		local d, h, d1, d2, d3, d1h1, d1h2, d2h1, h1, h2, h3 = unpack(coefients)
		d, h = d-x, h-y
		return d1*d + d2*d^2 +d3*d^3 + d1h1*d*h + d1h2*d*h^2 + d2h1*d^2*h + h1*h + h2*h^2 + h3*h^3
	end


	function find_bounding_triangle(u, v)
		local ix, iy = initialtangent(u[1], v[1], u[2], v[2], u[3], v[3], u[4], v[4])
		local fx, fy = finaltangent(u[1], v[1], u[2], v[2], u[3], v[3], u[4], v[4])
			local m1, m2 = (iy/ix), (fy/fx)
		if (ix == 0) then
			px = u[1]
			py = m2*px - m2*u[4] + v[4]
		elseif (fx == 0) then
			px = u[4]
			py = m1*px - m1*u[1] + v[1]
		else
			px = (m2*u[4] - m1*u[1] - v[4] + v[1])/(m2-m1)
			py = m1*px - m1*u[1] + v[1]
		end
		--print ("triangle", u[1], px, u[4], v[1], py, v[4])
		return {["x"] = {u[1], px, u[4]}, ["y"] = {v[1], py, v[4]}}
	 end


	function winding_implicit_quadratic(e, x, y)
		if (x < e.xmin) then return e.sign end
		local t = e.triangle
		sign = e.signx*e.sign
		if inside_triangle_implicit(t, x, y) then
			local px, py, pw = e.klm:apply(x, y, 1)
			if e.triangletoright then 
				implicit_test = (px^2 - py*pw) < 0
			else
				implicit_test = (px^2 - py*pw) >= 0
			end
			if (implicit_test) then
				return e.sign
			end
		else 
			local a, b, c = e.a[1], e.a[2], e.a[3]
			if e.triangletoright then 
				lineartest =  a*x + b*y + c <= 0
			else
				lineartest =  a*x + b*y + c < 0
			end
			if lineartest then
				return e.sign
			end
		end
		return 0
	end

	function winding_implicit_cubic(e, x, y)
		if (x < e.xmin) then return e.sign end
		local t = e.triangle
		if inside_triangle_implicit(t, x, y) then
			local coef = e.resultant_coefficients
			local res = implicit_resultant(coef, x, y)
			local curve_sign = e.intersectionsign
			local resultant_sign = util.sign(res)
			if e.triangletoright then 
				implicit_test = curve_sign == resultant_sign
			else
				implicit_test = curve_sign ~= resultant_sign
			end
			if (implicit_test) then
				return e.sign
			end
		else
			local a, b, c = e.a[1], e.a[2], e.a[3]
			if e.triangletoright then 
				lineartest =  a*x + b*y + c <= 0
			else
				lineartest =  a*x + b*y + c < 0
			end
			if lineartest then
				return e.sign
			end
		end
		return 0
	end

	function inside_triangle_implicit(t, px, py)
		-- using the baricentric test
		local d, na, nb, s
		local x, y = t.x, t.y
		x21, x31, qx = (x[2] - x[1]), (x[3] - x[1]), (px - x[1])
		y21, y31, qy = (y[2] - y[1]), (y[3] - y[1]), (py - y[1])
		d  = det(x21, x31, y21, y31)
		na = det(qx, x31, qy, y31)
		nb = det(x21, qx, y21, qy)
		s = util.sign(d)
		return (s*na > 0 and s*nb > 0 and s*(na + nb) <= s*d)
	end

	function inside_path(p, x, y, fill_type)
		local result, sign = 0, 0
		for i, e in pairs(p.segments) do
	 		if e.ymin < y and y <= e.ymax and e.xmax >= x then 
				sign = e.foo(e, x, y, i)
				result = result + sign
			end
		end
		--if fill_type == "eofill" then
		--	return (math.abs(result) % 2 ~= 0)
		--elseif fill_type == "fill" then 
			return result ~= 0
		--end
		--return true
	end

-- Sampling
	--used before blue
	function semi_random_grid(size)
		local init = -0.5 
		local dxy = 1/size
			local position = {}
		for j = 1, size do
			position[j] = {}
			for i = 1, size do
				local rx, ry = math.random() - 0.5, math.random() - 0.5
				px, py = init + dxy*(i - 0.5 + rx) , init + dxy*(j - 0.5 + ry)
				position[i + size*(j-1)] = {px, py}
				last = i + size*(j-1)
			end
		end
		return position
	end

	function supersampling(colors, n)
		local color = ungamma(colors[1])
		for i = 2, n do 
			color = adds(color, ungamma(colors[i]))
		end
		return gamma(mult(color, 1/n))
	end
	
	local g = 2.2
	
	function gamma(color)

		return {color[1]^g, color[2]^g, color[3]^g, color[4]^g} 
	end

	function ungamma(color)
		return {color[1]^(1/g), color[2]^(1/g), color[3]^(1/g), color[4]^(1/g)} 
	end

-- Basics
	local function preparescene(scene)
		scene.inverse_shapes_list = {}
		for i, element in ipairs(scene.elements) do
			if element.shape.type == "circle" then
				prepare_circle(scene, element.shape)
			elseif element.shape.type == "triangle" then	
				prepare_triangle(scene, element.shape)
			elseif element.shape.type == "polygon" then
				prepare_polygon(scene, element.shape)		
			elseif element.shape.type == "path" then
				prepare_path(scene, element.shape)
			end
			if element.paint.type == "lineargradient" then
				prepare_gradient_linear(scene, element.paint)
			elseif element.paint.type == "radialgradient" then
				prepare_gradient_radial(scene, element.paint)
			elseif element.paint.type == "solid" then
				element.paint.data[4] = element.paint.data[4]*element.paint.opacity
				element.paint.foo = coloring_solid
			elseif element.paint.type == "texture" then
				prepare_texture(scene, element.paint, element.shape)
				element.paint.foo = coloring_texture
			end
			scene.inverse_shapes_list[scene.size - i + 1] = element		
		end	
			
		-- GRID TEST --
		local svgout = require("export_cell")
		GRID = prepareGrid(scene)
		
		--[[
		for i = 1,#cells do
			for j = 1,#cells[i] do
				svgout.export_cell(cells[i][j], "cell" .. i .. j)
			end
		end ]]

		return scene
	end

	-- verifies that there is nothing unsupported in the scene
	local function checkscene(scene)
	    for i, element in ipairs(scene.elements) do
	        assert(element.type == "fill" or
	        	element.type == "eofill", "unsupported element")
	        assert(element.shape.type == "circle" or
	               element.shape.type == "triangle" or
		    	element.shape.type == "path" or
	               element.shape.type == "polygon", "unsuported primitive")
	        assert(element.paint.type == "solid" or
	        	element.paint.type == "lineargradient" or
	        	element.paint.type == "radialgradient" or
	        	element.paint.type == "texture", "unsupported paint")
			scene.size = i
	    end
	end

	-- output formatted string to stderr
	local function stderr(...)
	    io.stderr:write(string.format(...))
	end

	function _M.render(scene, viewport, file)	
		--dump(scene)
		print("Default: supersampling option: 1, gamma 2.2")
		option = 1
		assert( option == 1 or option == 8 or option == 16 or option == 32 or option == 64, 
			"Option not supported. Options: 1, 8, 16, 32, 64")
		local time = chronos.chronos()
		    -- make sure scene does not contain any unsuported content
		    checkscene(scene)
		    -- transform and prepare scene for rendering
		    local vxmin, vymin, vxmax, vymax = unpack(viewport, 1, 4)
		    local width, height = vxmax-vxmin, vymax-vymin
		    scene.width, scene.height =  width, height
		    scene = preparescene(scene)
		    -- get viewport
		stderr("preprocess in %.3fs\n", time:elapsed())
		time:reset()
		    -- get image width and height from viewport
		    -- allocate output image
		    local img = image.image(width, height)
		    -- render
		    for i = 1, height do
		stderr("\r%5g%%", floor(1000*i/height)/10)
		        local y = vymin+i-1.+.5
		        for j = 1, width do
		            local x = vxmin+j-1.+.5
		            if option > 1 then
		            	local newcolors, dp = {}, blue[option]
			            for k = 1, option*2, 2 do
			            	newcolors[(k+1)/2] = sample(scene, x+dp[k], y+dp[k+1])
			            end
			           	color = supersampling(newcolors,option)
			        else 
			        	color = sample(scene, x, y)
			        end--]]
		            img:set(j, i, unpack(color, 1, 4))
		        end
		    end
		stderr("\n")
		stderr("rendering in %.3fs\n", time:elapsed())
		time:reset()
		    -- store output image
		    image.png.store8(file, img)
		stderr("saved in %.3fs\n", time:elapsed())
	end

	function det(x11, x12, x21, x22)
		return (x11*x22 - x12*x21) 
	end

	function atan2(dy, dx)
		if (dy < 0) then
			if (dx < 0) then return math.pi - math.atan(math.abs(dy/dx))
			elseif (dx > 0) then return math.atan(math.abs(dy/dx))
			else return math.pi/2 end
		elseif (dx < 0) then return math.pi - math.atan(math.abs(dy/dx))
		elseif (dx > 0) then return math.atan(math.abs(dy/dx))
		end
		return math.pi/2
	end

	function mult(tb, cons)
		local table = {}
		for i, e in ipairs(tb) do
			table[i] = cons*e
		end
		return table
	end


	function adds(tb1, tb2)
		local table = {}
		for i, e in ipairs(tb1) do
			table[i] = tb1[i] + tb2[i]
		end
		return table
	end
	return _M
