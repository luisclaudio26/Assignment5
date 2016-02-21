local driver = require"driver"
local image = require"image"
local chronos = require"chronos"

local unpack = table.unpack
local bezier = require("bezier")
local floor, ceil, min, max, sqrt, inf = math.floor, math.ceil, math.min, math.max, math.sqrt, math.huge

local _M = driver.new()
local prepare = {}

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
		local x, y =  event_list[#event_list][2], event_list[#event_list][3]
		for i = #event_list -1, 1, -1 do
			next_y = event_list[i][3]
			if next_y == y then
				next_x = event_list[i][2]
			else
				next_x = 1
			end
			if w ~= 0 then
				for dx = next_x, x-1 do
					cells[dx][y].initialWindingNumber = w
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
			x, y = event_list[i][3], next_y
		end
		if w ~= 0 then
			for dx = 1, x do
				cells[dx][y].initialWindingNumber = w
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
    			["initialWindingNumber"] = 0, ["segments"] = {}, ["border"] = false }

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
		return {["x"] = x, ["y"] = y, ["foo"] = winding_number_linear,
			["sign"] = -1}
	end
	-----------------------------------------------------------------------------------

	--DONE AND WORKING
	local function walkInPath(path, cells)
	    -- For each segment inside the path, walk through it in a Brenseham/Tripod-fashion
	    -- push segment into intersecting cell
	    -- write events to event_list
	    for l, segment in ipairs(path.segments) do
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
		    	table.insert(cells[i][j].segments, segment)
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
		    			table.insert(event_list, {1, i, j})
	    				--print(i, j, "->", i, j+1)
		    			i, j = i, j+1
		    		end
		    	elseif intersectSegmentCell(i, j-1, cells, segment) then -- leave under, decrement initial winding number
	    			flag = true
		    	--	print("casedown")
	    			table.insert(event_list, {-1, i, j})
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
			    				local s = createClosingCellSegment(path, {x[#x], cells[i][j].ymax}, {x[#x], y[#y]})
			    				table.insert(cells[i][j].segments, s)
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
				walkInPath(e.shape, cells)
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
	end


	--TODO
	local function export_cell(i,j,cell)
	    -- 4) compose scene
	    -- 5) Export svg

	    -- RETURN: VOID, but exports a .svg file
	end



---------------------- shortcut operations--------------------------------------
local function sign(x) return (x<0 and -1) or 1 end

function table.empty (self)
    for _, _ in pairs(self) do  return false end
    return true
end

local function truncate_parameter(t)
    if t < 0 or t == -math.huge then t = 0 --????
    elseif t > 1 or t == math.huge then t = 1
    end
    return t
end

-- I use this to transform and add to the bounding box m points
local function dataTransform(shape,data,pos,m)
  for i=0,(m-1) do
    data[pos+2*i], data[pos+(2*i+1)] = shape.xf : apply(data[pos+2*i],data[pos+(2*i+1)])
    shape.xmin = min(shape.xmin or  inf, data[pos+2*i])
  	shape.ymin = min(shape.ymin or  inf, data[pos+2*i+1])
  	shape.xmax = max(shape.xmax or -inf, data[pos+2*i])
  	shape.ymax = max(shape.ymax or -inf, data[pos+2*i+1])
  end
  return data
end

---------------------------Auxiliarie functions---------------------------------
local function blend(f,b)
  return {f[4]*f[1]+(1-f[4])*b[4]*b[1],f[4]*f[2]+(1-f[4])*b[4]*b[2],f[4]*f[3]+(1-f[4])*b[4]*b[3],f[4]+(1-f[4])*b[4]}
end

local function cross_with_e1(v1, v2, v3)
  local a, b, c, d = v1[1], v1[2], v1[3], v1[4]
  local l, e, f, g = v2[1], v2[2], v2[3], v2[4]
  return -d*f + c*g,d*e - b*g,-c*e + b*f
end

local function bezier_to_poly3(x0,y0,x1,y1,x2,y2,x3,y3)
  return {x0, -3*x0 + 3*x1, 3*x0 - 6*x1 + 3*x2, -x0 + 3*x1 - 3*x2 + x3},
         {y0, -3*y0 + 3*y1, 3*y0 - 6*y1 + 3*y2, -y0 + 3*y1 - 3*y2 + y3}
end

local function difP(p)
  local dp = {0}
  for i=1,(#p-1) do dp[i] = i*p[i+1] end
  return dp
end

local function multP(p,q)
  local n, m, h = #p, #q, {}
  for i=1,(n+m-1) do h[i] = 0 end
  for i=1,n do
    for j=1,m do
      h[i+j-1] = h[i+j-1] + p[i]*q[j]
    end
  end
  return h
end

local function compute_rational_maxima(x0, x1, x2, w)
    local a = 2*(-1 + w)*(x0 - x2)
    local b = 2*(x0 - 2*w*x0 + 2*x1 - x2)
    local c = 2*(w*x0 - x1)

    n, r1, s1, r2, s2 = _M.quadratic(a, b, c)

    local out1, out2 = 0, 0
    if n > 0 then out1 = r1/s1 end
    if n > 1 then out2 = r2/s2 end

    out1, out2 = truncate_parameter(out1), truncate_parameter(out2)
    return out1, out2
end

local bezier_cut = {0,bezier.cut2,bezier.cut3}

local function polyToFun(p)
  return function (t)
    local n = (#p-1)
    local v = p[1]
    for i=1,n do v = v + p[i+1]*(t^i) end
    return v
  end
end

local function newton_aux(f,df,t0,t1,t,dt,y,dy,it)
  if math.abs(dt) < 10^(-12) or it > 50 or y == 0 then
    return t
  else
    if type(dy) ~= "number" or dy == math.huge or ((t-t1)*dy-y)*((t-t0)*dy-y) > 0 or math.abs(2*y) > math.abs(dy*dt) then
      -- bisection
      tm = (t0+t1)*(1/2); dtm = (t1-t0)*(1/2)
      ym = f(tm); dym = df(tm)
      if ym > 0 then
        return newton_aux(f,df,t0,tm,tm,dtm,ym,dym,it+1)
      else
        return newton_aux(f,df,tm,t1,tm,dtm,ym,dym,it+1)
      end
    else
      -- newton step
      dtm = -y*(1/dy); tm = t + dtm
      ym = f(tm); dym = df(tm)
      if ym > 0 then
        return newton_aux(f,df,t0,tm,tm,dtm,ym,dym,it+1)
      else
        return newton_aux(f,df,tm,t1,tm,dtm,ym,dym,it+1)
      end
    end
  end
end

local function safe_newton_raphson(f,df,a,b)
  local y0, y1, x = f(a), f(b), (a+b)*.5
  assert(y0*y1 <= 0, "there is no confirmation of root on the interval")
  if y0 > 0 then
    return newton_aux(f,df,b,a,x,b-a,f(x),df(x),1)
  else
    return newton_aux(f,df,a,b,x,b-a,f(x),df(x),1)
  end
end

local function poly_num_root(p,a,b)
  local Dp = difP(p)
  local funP, fundP = polyToFun(p), polyToFun(Dp)
  local So, S = roots[math.min(#Dp-1,4)](Dp,a,b), {}
  for i,t in pairs(So) do
    if t>a and t<b then S[#S+1] = t end
  end
  local roots = {}
  if funP(a) == 0 then roots[#roots+1] = a end
  if funP(b) == 0 then roots[#roots+1] = b end
  if #S==0 then
    if funP(a)*funP(b) > 0 then return {} end
    return {safe_newton_raphson(funP,fundP,a,b)}
  end
  if funP(a)*funP(S[1]) < 0 then roots[1] = safe_newton_raphson(funP,fundP,a,S[1])
  elseif funP(S[1])==0 then roots[#roots+1] = S[1]
  end
  for i=1,(#S-1) do
    if funP(S[i])*funP(S[i+1]) < 0 then roots[#roots+1] = safe_newton_raphson(funP,fundP,S[i],S[i+1])
    elseif funP(S[i+1]) == 0 then roots[#roots+1] = S[i+1]
    end
  end
  if funP(S[#S])*funP(b) < 0 then roots[#roots+1] = safe_newton_raphson(funP,fundP,S[#S],b)
  elseif funP(S[#S])==0 then roots[#roots+1] = S[#S]
  end
  return roots
end

roots = {
  function (p) --degree 1
    if p[2] then return -p[1]/p[2] end
  end,
  function (p) --degree 2
    local n,t1,s1,t2,s2 = _M.quadratic(p[3],p[2],p[1])
    if n==2 then return {t1*(1/s1),t2*(1/s2)}
    elseif n==1 then return {t1*(1/s1)}
    else return {}
    end
  end,
  function (p) --degree 3
    local n,t1,s1,t2,s2,t3,s3 = _M.cubic(p[4],p[3],p[2],p[1])
    if n==3 then return {t1*(1/s1),t2*(1/s2),t3*(1/s3)}
    elseif n==2 then return {t1*(1/s1),t2*(1/s2)}
    elseif n==1 then return {t1*(1/s1)}
    else return {}
    end
  end,
  function (p,a,b) --any degree
    local a, b = a or -math.huge, b or math.huge
    return poly_num_root(p,a,b)
  end
}

local function biP_to_fun(p)
  local n = #p
  return function(x,y)
    local v = 0
    for i=1,n,3 do
      local dx, dy, a = p[i], p[i+1], p[i+2]
      v = v + a*(x^dx)*(y^dy)
    end
    return v
  end
end

local function inside(shape,x,y)
  return x >= shape.xmin and x <= shape.xmax and y >= shape.ymin and y <= shape.ymax
end

-------------------------------Winding Number-----------------------------------
local winding = {}

function winding.line(c,x,y)
  if y > c.ymin and y <= c.ymax and c.l(x,y) < 0 then
    return c.s
  else return 0 end
end


local function windingNumber (path,x,y)
  local wN = 0
  for i,segment in ipairs(path) do
    wN = wN + winding[segment.type](segment,x,y)
  end
  return wN
end



-------------------------------Implicitization----------------------------------
local implicitEq = {
  function (x0,y0,x1,y1)
    local a, b = y1-y0, x0-x1
    local c, s = -(a*x0 + b*y0), sign(a)
    local p = {0,0,s*c,0,1,s*b,1,0,s*a}
    return biP_to_fun(p), s
  end
}

local implicit = {
  function (x0,y0,x1,y1)
    local l,s = implicitEq[1](x0,y0,x1,y1)
    return {x = {x0,x1}, y = {y0,y1}, ymin = min(y0,y1), ymax = max(y0,y1), l = l, s = s, type = "line"}
  end
}


-- prepare paths
-------------------------------Instructions-------------------------------------
local execute = {}
function execute.begin_closed_contour(shape)
  shape.beg = shape.pos + 1
  dataTransform(shape,shape.data,shape.beg,1)
end

execute.begin_open_contour = execute.begin_closed_contour

function execute.degenerate_segment(shape) shape.degenerated = true end

function execute.end_closed_contour(shape)
  shape.segments[#shape.segments+1] = implicit[1](shape.data[shape.pos],shape.data[shape.pos+1],shape.data[shape.beg],shape.data[shape.beg+1])
end

function execute.end_open_contour(shape)
  if not shape.degenerated then
    shape.degenerated = false
  else
    shape.segments[#shape.segments+1] = implicit[1](shape.data[shape.pos],shape.data[shape.pos+1],shape.data[shape.beg],shape.data[shape.beg+1])
  end
end

function execute.linear_segment(shape)
  local data, pos = shape.data, shape.pos
  dataTransform(shape,data,pos+2,1)
  shape.segments[#shape.segments+1] = implicit[1](data[pos],data[pos+1],data[pos+2],data[pos+3])
end

function execute.quadratic_segment(shape)
  local data, pos = shape.data, shape.pos
  dataTransform(shape,data,pos+2,2)
  local x0, y0 = data[pos], data[pos+1]
  local x1, y1 = data[pos+2], data[pos+3]
  local x2, y2 = data[pos+4], data[pos+5]

  local partition, v = {}, 1./(x0 - 2*x1 + x2)
  partition[1], partition[4] = 0, 1
  partition[2] = truncate_parameter( (x0-x1) * v )
  partition[3] = truncate_parameter( (y0-y1) * v )
  table.sort( partition )

  for i = 2, 4 do
      if partition[i-1] ~= partition[i] then
          u0, v0, u1, v1, u2, v2 = bezier_cut[2](partition[i-1], partition[i], x0, y0, x1, y1, x2, y2)
          shape.segments[#shape.segments+1] = implicit[2](u0, v0, u1, v1, u2, v2)
      end
  end
end

function execute.cubic_segment(shape)
  local data, pos = shape.data, shape.pos
  dataTransform(shape,data,pos+2,3)

  --Finding inflections and double points
  local x0, y0 = data[pos], data[pos+1]
  local x1, y1 = data[pos+2], data[pos+3]
  local x2, y2 = data[pos+4], data[pos+5]
  local x3, y3 = data[pos+6], data[pos+7]

  local px, py = bezier_to_poly3(x0,y0,x1,y1,x2,y2,x3,y3)
  local d2,d3,d4 = cross_with_e1(px,py)

  local partition, degree = {0,1}, 3
  if d2==0 and d3==0 then
    if d4==0 then--Case it degenerates to a straight line
      shape.segments[#shape.segments+1] = implicit[1](x0,y0,x3,y3)
      return shape
    else --Case it degenerates to a quadratic
      degree = 2
      local xaux = (3/2)*(x1 - x0) + x0
			local yaux = (3/2)*(y1 - y0) + y0
      local w = 1./(x0 - 2*xaux + x2)
      partition[3] = truncate_parameter( (x0-xaux) * w )
      partition[4] = truncate_parameter( (y0-yaux) * w )
    end
  else
    local p = {-d4, 3*d3, -3*d2}
    local q = {d3^2 - d2*d4, - d2*d3, d2^2}
    local inflections, doublePts = roots[2](p), roots[2](q)
    local P = multP(difP(px),difP(py))
    local monotone_partition = roots[min(#P-1,4)](P,0,1)
    for i,e in ipairs(inflections) do partition[#partition+1] = truncate_parameter(e) end
    for i,e in ipairs(doublePts) do partition[#partition+1] = truncate_parameter(e) end
    for i,e in ipairs(monotone_partition) do partition[#partition+1] = truncate_parameter(e) end
  end
  table.sort(partition)

  for i = 2, #partition do
      if partition[i-1] ~= partition[i] then
        local u0, v0, u1, v1, u2, v2, u3, v3 = bezier_cut[degree](partition[i-1], partition[i], x0, y0, x1, y1, x2, y2, x3, y3)
        shape.segments[#shape.segments+1] = implicit[degree](u0, v0, u1, v1, u2, v2, u3, v3)
      end
  end
end

function execute.rational_quadratic_segment(shape)
  local data, pos = shape.data, shape.pos
  data[pos+2], data[pos+3], data[pos+4] = shape.xf : apply(data[pos+2], data[pos+3], data[pos+4])
  data[pos+5], data[pos+6] = shape.xf : apply(data[pos+5], data[pos+6])

  local x0, y0 = data[pos], data[pos+1]
  local x1, x1, w = data[pos+2], data[pos+3], data[pos+4]
  local x2, y2 = data[pos+5], data[pos+6]
  local invW = 1./w
  shape.xmin = min(shape.xmin, x0, x1*invW, x2)
  shape.ymin = min(shape.ymin, y0, y1*invW, y2)
  shape.xmax = max(shape.xmax, x0, x1*invW, x2)
  shape.ymax = max(shape.ymax, y0, y1*invW, y2)

  -- Find maxima
  t = {}
  t[1], t[6] = 0, 1
  t[2], t[3] = compute_rational_maxima(x0, x1, x2, w)
  t[4], t[5] = compute_rational_maxima(y0, y1, y2, w)
  table.sort( t )

  for i = 2, 6 do
      if t[i-1] ~= t[i] then
          local u0, v0, u1, v1, r, u2, v2 = bezier.cut2rc(t[i-1], t[i], x0, y0, x1, y1, w, x2, y2)
          shape.segments[#shape.segments+1] = implicit[d](u0, v0, u1, v1, u2, v2, r)
      end
  end
end

--------------------------------------------------------------------------------
function prepare.path(shape)
  shape.segments = {}
  shape.beg, shape.degenerated = 1, false
  for j, instruction in pairs(shape.instructions) do
    shape.pos = shape.offsets[j]
    execute[instruction](shape)
  end
end

-- prepare for simple paths
function prepare.triangle(shape)
  local x1,y1,x2,y2,x3,y3 = unpack(dataTransform(shape,{shape.x1,shape.y1,shape.x2,shape.y2,shape.x3,shape.y3},1,3))
  shape.segments = {implicit[1](x1,y1,x2,y2), implicit[1](x2,y2,x3,y3), implicit[1](x3,y3,x1,y1)}
  return shape
end

function prepare.polygon(shape)
  local D = shape.data
  local n = #D
  shape.segments = {}
  dataTransform(shape,D,1,n*.5)
  for i=1,(n-3),2 do
    shape.segments[#shape.segments+1] = implicit[1](D[i],D[i+1],D[i+2],D[i+3])
  end
  shape.segments[#shape.segments+1] = implicit[1](D[n-1],D[n],D[1],D[2])
end

function prepare.circle(shape)
  -- it is formed by 4 arcs covering each quarter of the unit circle (Using Four Vertex Theorem[Manfredo C.] )
  local s = sqrt(2)/2          -- sin(pi/2)
  local cx, cy, r = shape.cx, shape.cy, shape.r
  local arc1 = shape.xf * _M.xform(-r+cx,-r*s+cx,   cx,
                                      cy, r*s+cy, r+cy,
                                      1,      s,    1)
  local arc2 = shape.xf * _M.xform( cx, r*s+cx, r+cx,
                                    r+cy, r*s+cy,   cy,
                                       1,      s,    1)
  local arc3 = shape.xf * _M.xform( r+cx, r*s+cx,   cx,
                                      cy,-r*s+cy,-r+cy,
                                       1,      s,    1)
  local arc4 = shape.xf * _M.xform(   cx,-r*s+cx,-r+cx,
                                   -r+cy,-r*s+cy,   cy,
                                       1,      s,    1)
  local x = {arc1[1],arc1[2]/arc1[8],arc1[3],
             arc2[1],arc2[2]/arc2[8],arc2[3],
             arc1[1],arc3[2]/arc3[8],arc3[3],
             arc1[1],arc4[2]/arc4[8],arc4[3]}
  local y = {arc1[4],arc1[2]/arc1[8],arc1[3],
             arc2[4],arc2[2]/arc2[8],arc2[3],
             arc1[4],arc3[5]/arc3[8],arc3[3],
             arc1[4],arc4[5]/arc4[8],arc4[3]}
  shape.xmin = min(unpack(x))
  shape.xmax = max(unpack(x))
  shape.ymin = min(unpack(y))
  shape.ymax = max(unpack(y))

  shape.segments = {implicit[2](arc1),implicit[2](arc2),implicit[2](arc3),implicit[2](arc4)}
end
--------------------------------------------------------------------------------
--------------------------Preparing Paint---------------------------------------
function prepare.solid(paint)
  paint.data[4] = paint.data[4]*paint.opacity
  paint.color = function (self,x,y)
    return self.data
  end
end






--------------------------------------------------------------------------------
-- prepare scene for sampling and return modified scene
local function preparescene(scene)
  for i, element in ipairs(scene.elements) do
    element.shape.xf = scene.xf * element.shape.xf
    prepare[element.shape.type](element.shape)
    element.paint.xf = element.paint.xf:inverse() * scene.xf:inverse()
    prepare[element.paint.type](element.paint)
  end
  -- prepare.grid(scene)
  return scene
end

-- sample scene at x,y and return r,g,b,a
local function sample(scene, x, y)
  local color = {1,1,1,1}
  for ind,e in ipairs(scene.elements) do
    local shape, paint = e.shape, e.paint
    if inside(shape,x,y) then
      local wN = windingNumber(shape.segments,x,y)
      if e.type == "fill" and wN ~= 0 then
        color = blend(paint:color(x,y),color)
      elseif e.type == "eofill" and wN % 2 == 1 then
        color = blend(paint:color(x,y),color)
      end
    end
  end
  return unpack(color,1,4)
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
    end
end

-- output formatted string to stderr
local function stderr(...)
    io.stderr:write(string.format(...))
end

function _M.render(scene, viewport, file)
local time = chronos.chronos()
    -- make sure scene does not contain any unsuported content
    checkscene(scene)
    -- get viewport
    local vxmin, vymin, vxmax, vymax = unpack(viewport, 1, 4)
    -- get image width and height from viewport
    local width, height = vxmax-vxmin, vymax-vymin
    -- transform and prepare scene for rendering
    scene.width, scene.height = width, height
    scene = preparescene(scene)

stderr("preprocess in %.3fs\n", time:elapsed())
time:reset()
    -- allocate output image
    local img = image.image(width, height)
    -- render2D/Lua/ubuntu-1.02$ ./lua process.lua ../../assign-4-1.0/render.lua ../../rvg-1.03/qua.rvg out.png

    for i = 1, height do
stderr("\r%5g%%", floor(1000*i/height)/10)
        local y = vymin+i-1.+.5
        for j = 1, width do
            local x = vxmin+j-1.+.5
            img:set(j, i, sample(scene, x, y))
        end
    end
stderr("\n")
stderr("rendering in %.3fs\n", time:elapsed())
time:reset()
    -- store output image
    image.png.store8(file, img)
stderr("saved in %.3fs\n", time:elapsed())
end

return _M
