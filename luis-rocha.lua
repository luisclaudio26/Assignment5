local driver = require"driver"
local image = require"image"
local chronos = require"chronos"

local bezier, quadratic = require("lua.bezier"), require("lua.quadratic")
local bernstein = require("lua.bernstein")
local xform, noise = require("xform"), require("blue")
local unpack, pack = table.unpack, table.pack
local max, min = math.max, math.min 
local floor, ceil = math.floor, math.ceil
local abs = math.abs

local _M = driver.new()

local BGColor = require("lua.color").rgb(1,1,1,1)
local epsilon, max_iteration, gamma_factor = 0.0000000001, 50, 2.2

local n_cells_x, n_cells_y = 10, 10

-----------------------------------------------------------------------------------------
-------------------------------- AUXILIAR FUNCTIONS -------------------------------------
-----------------------------------------------------------------------------------------
local function sign(v)
    if v < 0 then return -1
    elseif v > 0 then return 1
    else return 0 end
end

local function transform_point(x, y, xf)
    local _x, _y, w = xf : apply(x, y, 1)
    return _x / w, _y / w
end

local function truncate_parameter(t)
    if t < 0 or t == -math.huge then t = 0 --????
    elseif t > 1 or t == math.huge then t = 1
    end
    return t
end

local function root_bisection(t0, t1, func, n_it)
    local n_it = n_it or 0
    local tm = (t0 + t1)*0.5

    -- Halting criterium
    local delta = t1 - t0
    if delta < epsilon then return tm end
    if n_it >= max_iteration then return false end

    -- Recursively subdivide
    local y = func(tm)
    local sy = sign(y)

    if sy == sign( func(t0) ) then
        return root_bisection(tm, t1, func, n_it + 1)
    elseif sy == sign( func(t1) ) then 
        return root_bisection(t0, tm, func, n_it + 1)
    elseif y == 0 then 
        return tm
    end
end

local function compute_cubic_maxima(x0, x1, x2, x3)
    local c = 3*(x1 - x0)
    local b = 6*(x0 - 2*x1 + x2)
    local a = 3*(-x0 + 3*x1 - 3*x2 + x3)

    local n, t1, s1, t2, s2 = quadratic.quadratic(a,b,c)
    local out1, out2 = 0, 0

    if n > 0 then out1 = t1/s1 end
    if n > 1 then out2 = t2/s2 end

    out1, out2 = truncate_parameter( out1 ), truncate_parameter( out2 )
    return out1, out2
end

local function compute_cubic_inflections(x0, y0, x1, y1, x2, y2, x3, y3, d2, d3, d4)
    local a, b, c = -3*d2, 3*d3, -d4

    local n, t1, s1, t2, s2 = quadratic.quadratic(a, b, c)
    local out1, out2 = 0, 0

    -- Other inflection is in infinity (which is consistent). Should it be
    -- truncated to 1?
    if n > 0 and s1 ~= 0 then out1 = t1/s1 end
    if n > 1 and s2 ~= 0 then out2 = t2/s2 end

    out1, out2 = truncate_parameter(out1), truncate_parameter(out2)

    return out1, out2
end

local function compute_cubic_doublepoint(x0, y0, x1, y1, x2, y2, x3, y3, d2, d3, d4)
    local a, b, c =  d2^2, -d2*d3, d3^2 - d2*d4
    local n, t1, s1, t2, s2 = quadratic.quadratic(a, b, c)
    local out1, out2 = 0, 0

    if n > 0 and s1 ~= 0 then out1 = t1/s1 end
    if n > 1 and s2 ~= 0 then out2 = t2/s2 end

    out1, out2 = truncate_parameter(out1), truncate_parameter(out2)

    return out1, out2
end

local function compute_rational_maxima(x0, x1, x2, w)
    local a = 2*(-1 + w)*(x0 - x2)
    local b = 2*(x0 - 2*w*x0 + 2*x1 - x2)
    local c = 2*(w*x0 - x1)

    n, r1, s1, r2, s2 = quadratic.quadratic(a, b, c)

    local out1, out2 = 0, 0
    if n > 0 then out1 = r1/s1 end
    if n > 1 then out2 = r2/s2 end

    out1, out2 = truncate_parameter(out1), truncate_parameter(out2)
    return out1, out2
end

local function alpha_composite(c1, a1, c2, a2)
    -- Assumes non-premultiplied values
    local num = c1*a1 + c2*a2*(1-a1)
    local den = a1 + a2*(1-a1)

    return num/den
end

local function search_in_ramp(ramp, value)
    -- Just a linear search. Other versions may sample
    -- the ramp and just do a look-up table for this
    for i = 1, #ramp-2, 2 do
        if ramp[i+2] >= value then return i end
    end
end

local function interpolate_colors(color1, color2, t)
    local out = {}
    for i = 1, 4 do
        out[i] = (1-t)*color1[i] + t*color2[i]
    end
    return out
end

local function fix_ramp(ramp)
    -- If 0.0 and 1.0 are not defined in the ramp, define it
    if ramp[1] ~= 0 then
        table.insert(ramp, 1, ramp[2]) -- Insert color
        table.insert(ramp, 1, 0) -- Insert offset
    end

    if ramp[#ramp-1] ~= 1 then
        local v = ramp[#ramp]
        table.insert(ramp, 1) -- Insert offset
        table.insert(ramp, v) -- Insert value
    end
end


local function compute_tangent_intersection(u0, v0, u1, v1, u2, v2, u3, v3, diagonal)
    -- We assume the curve to be MONOTONIC

    if u0 > u3 then
        -- Thanks Lua for making this possible
        u0, v0, u3, v3 = u3, v3, u0, v0
        u1, v1, u2, v2 = u2, v2, u1, v1
    end

    local outx, outy

    if u1 == u2 and v1 == v2 then
        -- First case: control points are coincident. Just return it
        return u1, v1
    elseif u0 == u1 and v0 == v1 then
        -- Second case: first point coincide with the second. Intersection
        -- will be then between line (u2,v2) -> (u3,v3) and x/y axis (if p2 is
        -- to the right/left of the diagonal linking p0 -> p3)
        local diag = diagonal(u2,v2)
        outx = u2
        outy = v2
    elseif u2 == u3 and v2 == v3 then
        -- Third case: dual to the second
        local diag = diagonal(u1,v1)
        outx = u1
        outy = v1
    else
        -- Fourth case: All points are different - intersection of tangents
        outx = (u0*(u3*(-v1 + v2) + u2*(v1 - v3)) + u1*(u3*(v0 - v2) + u2*(-v0 + v3)))
        outx = outx / (-(u2 - u3)*(v0 - v1) + (u0 - u1)*(v2 - v3))

        outy = (u3*(v0 - v1)*v2 + u0*v1*v2 - u2*v0*v3 - u0*v1*v3 + u2*v1*v3 + u1*v0*(-v2 + v3))
        outy = outy / (-(u2 - u3)*(v0 - v1) + (u0 - u1)*(v2 - v3))
    end

    return outx, outy
end

local function compute_implicit_line(x0, y0, x1, y1)
    local a = y1 - y0
    local b = x0 - x1
    local c = - a * x0 - b * y0
    local dysign = sign(a)

    return a*dysign, b*dysign, c*dysign
end

local function gamma(color)
    local out = {}
    for k,v in ipairs(color) do
        out[k] = v^gamma_factor
    end
    return out
end

local function ungamma(color)
    local f = 1/gamma_factor
    local out = {}
    for k,v in ipairs(color) do
        out[k] = v^f
    end
    return out
end

-------------------------------------------------------------------------------------
-------------------------------- GRID FUNCTIONS -------------------------------------
-------------------------------------------------------------------------------------
-- Contains (1) the cell coordinates and (2) the initial winding number increment
local event_list = {}

local function fixLineWindingNumber(line, start_i, end_i)
    -- Loop through line changing winding number 'till the net winding number is zero
    -- RETURN VOID
end

local function computeGridDimension(scene)
    -- this should return a "optimal" width and height after
    return n_cells_width, n_cells_height
end

local function makeGrid(window_width, window_height, n_cells_x, n_cells_y)
    -- returns an empty grid (a bidimensional table), where each cell 
    -- contains (1) its bounding box and (2) a table with the intersecting segments
    -- and the (3) initial winding number

    -- RETURN: A TABLE WITH FORMAT CELL[i][j] = {xmin, ymin, xmax, ymax, initialWindingNumber, segments = {} }

    local cell_w, cell_h = window_width/n_cells_x, window_height/n_cells_y
    local grid = {}

    grid.width, grid.height = n_cells_x, n_cells_y
    grid.cell_w, grid.cell_h = cell_w, cell_h
    grid.cells = {}

    for i = 1, n_cells_x do        
        grid.cells[i] = {}
        
        for j = 1, n_cells_y do
            local cell = grid.cells[i][j]

            -- CONVENTION: assumes box is closed in the left/bottom side, 
            -- open in the top/right side
            cell.xmin, cell.ymin = i*cell_w, j*cell_h
            cell.xmax, cell.ymax = (i+1)*cell_w, (j+1)*cell_h

            cell.initialWindingNumber = 0

            -- !!! Shapes must contain triplets in the form {segment = , fill_type = , paint = } !!!
            -- Or even maybe fill_type(segment, paint), just like a normal Element
            cell.shapes = {}
        end
    end

    return grid
end

local function intersectSegmentCell(x0, y0, x1, y1, segment)
    -- 1) Check path bounding box against Cell
    --      -> if it is fully inside, then return true
    --      -> if it is not, check if one of the extreme control points are inside the cell.
    -- 2) Ray cast -> check for intersection with right side
    --      -> If it does not intersect, check for intersection with top side (by rotating the cell)

    -- RETURN : BOOLEAN
    -- x(t) = xmin, x(t) = xmax, y(t) = ymin, y(t) = ymax 

    -- TODO: E quando o segmento cruzar dois lados da bounding box?
    -- A ordem dessas checagens deve estar correta!!!

    local direction

    -- Check boundaries
    local f = root_bisection(0, 1, function(t) return (segment.atx(t) - xmax) end )
    if f ~= false then direction = "right" goto orientation end

    f = root_bisection(0, 1, function(t) return (segment.atx(t) - xmin) end )
    if f ~= false then direction = "left" goto orientation end

    f = root_bisection(0, 1, function(t) return (segment.aty(t) - ymin) end )
    if f ~= false then direction = "bottom" goto orientation end

    f = root_bisection(0, 1, function(t) return (segment.aty(t) - ymax) end )
    if f ~= false then direction = "top" goto orientation end

    -- If we reach this part, segment is trapped inside box
    return "none"

    ::orientation::
    if segment.x0 >= x0 and segment.x0 < x1 and segment.y0 >= y0 and segment.y0 < y1 then
        return "leaving", direction
    else
        return "entering", direction
    end
end

local function walkInPath(element, grid)
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

    local shape, prim = element.shape, element.shape.primitives
    
    local seg_i = 1
    local i, j = getCell(prim[seg_i].x0, prim[seg_i].y0)
    local start_i, start_j = i, j

    repeat
        -- Get current cell and segment
        local cell = grid[i][j]
        local seg = prim[seg_i]

        -- Check whether segment is entering/leaving cell,
        -- and in what direction
        enterleave, direction = intersectSegmentCell(cell.xmin, cell.ymin, cell.xmax, cell.ymax, seg)

        -- Push segment to cell's list
        table.insert( cell.shapes, {segment = seg, fill_type = element.type, paint = element.paint} )

        local last_x, last_y = seg.last_point() -- Assumes untransformed
        if direction == "right" then
            j = j + 1

            --push straight line
            if last_y < cell.ymax then

                local line = {}
                if enterleave == "leaving" then
                    line.x0, line.y0, line.x1, line.y1 = last_x,last_y,last_x,cell.ymax
                else
                    line.x0, line.y0, line.x1, line.y1 = last_x,cell.ymax,last_x,last_y
                end

                local a, b, c = compute_implicit_line(line.x0, line.y0, line.x1, line.y1)
                line.implicit = function(x, y)
                    return a*x + b*y + c
                end

                table.insert(cell.shapes, {segment = line, fill_type = element.type, paint = element.paint})
            end
        elseif direction == "bottom" then
            i = i - 1

            if enterleave == "entering" then
                cell.initialWindingNumber = cell.initialWindingNumber + 1
            else
                cell.initialWindingNumber = cell.initialWindingNumber - 1
            end
        elseif direction == "left" then
            j = j - 1
        elseif direction == "top" then
            i = i + 1
        else 
            seg_i++ -- Segment is trapped inside cell
        end

    while not (i == start_i and j == start_j) and seg_i <= #prim

    -- TODO: Se o segment está entrando e acaba ali dentro, temos de pegar o próximo

    -- > Verifique se o segment sai pela direita, por cima, por baixo, pela esquerda ou se
    -- está totalmente dentro usando usando intersectSegmentCell
    -- > Insira o segment na lista de segmentos da célula
    -- > Dependendo do resultado, tome a próxima célula ([i+1, j], [i, j+1], dependendo do caso) e repita
    -- o processo ATÉ chegar na célula de partida
    -- >> Cuidado com o caso em que a célula sai/entra pela esquerda: neste caso temos que adicionar a reta
    -- extra.
    -- >> Também, a cada célula que sai/entra por baixo, devemos anotar na event_list
end

local function prepareGrid(rvg)
    -- 1) Compute grid dimensions
    -- 2) Create grid
    -- 3) loop through paths inside scene
    -- 4) Sort event_list -> insertion_sort (or any other stable sort)
    -- 5) fix initial winding numbers
    -- RETURN: FILLED GRID

    -- !!! THIS ASSUMES PREPROCESSED RVG !!!
    local window_w = rvg.viewport.xmax - rvg.viewport.xmin
    local window_h = rvg.viewport.ymax - rvg.viewport.ymin
    local grid = makeGrid(window_w, window_h, n_cells_x, n_cells_y)
    
    local scene = rvg.scene
    for i, element in ipairs(scene) do
        walkInPath(element, grid)
    end

    -- TODO:
    -- 4) Sort event_list -> insertion_sort (or any other stable sort)
    -- 5) fix initial winding numbers
end

local function getCell(x, y, grid)
    -- Compute coordinates of cell containing (x,y)
    -- RETURN: i, j -> INTEGERS!!!
    return ceil(x/grid.cell_w), ceil(y/grid.cell_h)
end

local function export_cell(cell)
    local paint, color = require("lua.paint"), require("lua.color")

    -- CREATE A TEST CELL --
    test_cell = {}
    test_cell.xmin, test_cell.xmax = 0, 100
    test_cell.ymin, test_cell.ymax = 0, 100
    test_cell.initialWindingNumber = 0
    
    test_cell.shapes = {{}, {}, {}, {}}

    test_cell.shapes[1].paint = paint.solid( color.rgb8(0,128,0) )
    test_cell.shapes[1].fill_type = "fill"
    test_cell.shapes[1].segment = {["type"] = "quadratic_segment", x0 = 0, y0 = 0, x1 = 50, y1 = 70, x2 = 100, y2 = 0}

    test_cell.shapes[2].paint = paint.solid( color.rgb8(0,0,255) )
    test_cell.shapes[2].fill_type = "fill"
    test_cell.shapes[2].segment = {["type"] = "linear_segment", x0 = 0, y0 = 0, x1 = 100, y1 = 0}

    test_cell.shapes[3].paint = paint.solid( color.rgb8(128,0,0) )
    test_cell.shapes[3].fill_type = "fill"
    test_cell.shapes[3].segment = {["type"] = "cubic_segment", x0 = 0, y0 = 0, x1 = 25, y1 = 80, x2 = 75, y2 = 80, x3 = 100, y3 = 0}    

    test_cell.shapes[4].paint = paint.solid( color.rgb8(0,100,100) )
    test_cell.shapes[4].fill_type = "fill"
    test_cell.shapes[4].segment = {["type"] = "rational_quadratic_segment", x0 = 0, y0 = 0, x1 = 50, y1 = 70, w = 2.0, x2 = 100, y2 = 0}

    -- This is how we call it. File will be output to a file called cell.svg
    require("export_cell").export_cell(test_cell)
end

-----------------------------------------------------------------------------------------
-------------------------------- PATH PREPROCESSING -------------------------------------
-----------------------------------------------------------------------------------------
local prepare_table = {}
prepare_table.instructions = {}
prepare_table.push_functions = {}

-------------------------------- TODO: (1) create atx(), aty() for all segments
------ Push primitives ---------       (2) create a function last_point() that returns the last control point of segment
--------------------------------            -> But points are transformed; what to do now?
function prepare_table.push_functions.linear_segment(x0, y0, x1, y1, holder, virtual)
    
    -- Virtual flag: a Virtual linear segment must be considered when filling
    -- the shape, but not if we're stroking. This is just a "fake" linear segment
    -- we use to close and open path when we need to fill it.
    
    local n = #holder + 1
    holder[n] = {}

    holder[n].type = "linear_segment"
    holder[n].virtual = virtual
    holder[n].x0, holder[n].y0 = x0, y0
    holder[n].x1, holder[n].y1 = x1, y1

    -- Precompute implicit equation and bounding box
    holder[n].a, holder[n].b, holder[n].c = compute_implicit_line(x0, y0, x1, y1)

    local xmin, xmax = min(x0, x1), max(x0, x1)
    local ymin, ymax = min(y0, y1), max(y0, y1)

    holder[n].dysign = sign( y1 - y0 )
    holder[n].xmin, holder[n].xmax = xmin, xmax
    holder[n].ymin, holder[n].ymax = ymin, ymax
end

function prepare_table.push_functions.degenerate_segment(x0, y0, dx0, dy0, dx1, dy1, holder)
    local n = #holder + 1
    holder[n] = {}

    holder[n].type = "degenerate_segment"
    holder[n].x0, holder[n].y0 = x0, y0
    holder[n].dx0, holder[n].dy0 = dx0, dy0
    holder[n].dx1, holder[n].dy1 = dx1, dy1
end

function prepare_table.push_functions.quadratic_segment(u0, v0, u1, v1, u2, v2, holder)

    -- Header info
    local n = #holder + 1 
    holder[n] = {}
    holder[n].type = "quadratic_segment"

    -- Bounding box (of untransformed points)
    local maxy, miny = max(v0, v2), min(v0, v2)
    local maxx, minx = max(u0, u2), min(u0, u2)
    holder[n].xmax, holder[n].xmin = maxx, minx
    holder[n].ymax, holder[n].ymin = maxy, miny

    -- Compute transformation to origin and translate control points
    local trans = xform.translate(-u0, -v0)

    u0, v0 = transform_point(u0, v0, trans)
    u1, v1 = transform_point(u1, v1, trans)
    u2, v2 = transform_point(u2, v2, trans)

    holder[n].dysign = sign(v2 - v0)

    -- Triangle test: compute the diagonal cutting
    -- the bounding box in two triangles
    local a, b, c = compute_implicit_line(u0, v0, u2, v2)

    holder[n].diagonal = function(x, y)
        return sign( a*x + b*y + c )
    end

    holder[n].mid_point_diagonal = holder[n].diagonal( bezier.at2(0.5, u0, v0, u1, v1, u2, v2) )

    -- Compute implicit equation based on resultant
    local det = xform.xform(u0, v0, 1, u1, v1, 1, u2, v2, 1) : det()
    local imp

    if det ~= 0 then
        local imp_sign = sign( 2*v2*(u1*v2 - u2*v1) )
        imp = function(x, y)
            local imp_sign = sign( 2*v2*(u1*v2 - u2*v1) )
            local diag1 = (-2*u1*y + 2*x*v1)*(2*u2*v1 - 2*u1*v2) 
            local diag2 = ((2*u1 - u2)*y + x*(-2*v1 + v2))^2

            local eval = diag2 - diag1
            if imp_sign < 0 then eval = -eval end

            return eval 
        end
    else
        imp = holder[n].diagonal
    end

    holder[n].implicit = imp

    -- Store transformed control points
    holder[n].x0, holder[n].y0 = u0, v0
    holder[n].x1, holder[n].y1 = u1, v1
    holder[n].x2, holder[n].y2 = u2, v2

    holder[n].scene_to_canonic = trans

    -- This will be useful
    holder[n].atx = function(t)
        return bernstein.lerp2(t, t, u0, u1, u2)
    end

    holder[n].aty = function(t)
        return bernstein.lerp2(t, t, v0, v1, v2)
    end

    holder[n].last_point = function() return u2, v2 end
end

function prepare_table.push_functions.cubic_segment(u0, v0, u1, v1, u2, v2, u3, v3, holder)
    
    -- Header
    local n = #holder + 1
    holder[n] = {}
    holder[n].type = "cubic_segment"

    -- Untransformed bounding box
    holder[n].xmax, holder[n].xmin = max(u0, u3), min(u0, u3)
    holder[n].ymax, holder[n].ymin = max(v0, v3), min(v0, v3)

    -- Translate first control point to the origin
    local trans = xform.translate(-u0, -v0)

    u0, v0 = transform_point(u0, v0, trans)
    u1, v1 = transform_point(u1, v1, trans)
    u2, v2 = transform_point(u2, v2, trans)
    u3, v3 = transform_point(u3, v3, trans)

    holder[n].scene_to_canonic = trans

    -- Triangle test: compute the diagonal cutting
    -- the bounding box in two triangles
    local a, b, c = compute_implicit_line(u0, v0, u3, v3)

    holder[n].diagonal = function(x, y)
        return sign( a*x + b*y + c )
    end

    holder[n].mid_point_diagonal = holder[n].diagonal( bezier.at3(0.5, u0, v0, u1, v1, u2, v2, u3, v3) )

    -- Triangle test: compute intersection of tangents
    local inter_x, inter_y = compute_tangent_intersection(u0, v0, u1, v1, u2, v2, u3, v3, holder[n].diagonal)

    local det22 = function(a, b, c, d)
        return a*d - b*c
    end

    local triangle_det = det22(inter_x-u0, u3-u0, inter_y-v0, v3 - v0)

    holder[n].inside_triangle = function(x, y)
        local a = det22(x - u0, u3 - u0, y - v0, v3 - v0)
        local b = det22(inter_x - u0, x - u0, inter_y - v0, y - v0)
        local d_sign = sign(triangle_det)

        return d_sign * a > 0 and 
            d_sign * b > 0 and 
            d_sign * (a + b) < d_sign * triangle_det
    end

    -- Compute implicit function (missing sign test, degenerate test)
    local det1 = xform.xform(u0,v0,1,u1,v1,1,u2,v2,1) : det()
    local det2 = xform.xform(u1,v1,1,u2,v2,1,u3,v3,1) : det()

    local imp

    if det1 ~= 0 or det2 ~= 0 then

        local a1 = (v1 - v2 - v3)
        local a2 = -(4*v1^2 - 2*v1*v2 + v2^2)*u3^2
        local a3 = (9*v2^2 - 6*v2*v3 - 4*v3^2)*u1^2
        local a4 = (9*v1^2 - 12*v1*v3 - v3^2)*u2^2
        local a5 = 2*u1*u3*(-v2*(6*v2 + v3) + v1*(3*v2 + 4*v3))
        local a6 = 2*u2*(u3*(3*v1^2 - v2*v3 + v1*(-6*v2 + v3)) + u1*(v1*(9*v2 - 3*v3) - v3*(6*v2 + v3)))
        local imp_sign = sign( a1*(a2+a3+a4+a5-a6) )

        imp = function(x, y)
            f1 = -(-9*u2*v1 + 3*u3*v1 + 9*u1*v2 - 3*u1*v3)
            f2 = (-3*u1*y + 3*x*v1)
            f3 = (-9*u2*v1 + 3*u3*v1 + 9*u1*v2 - 3*u1*v3)
            f4 = ((6*u1 - 3*u2)*y + x*(-6*v1 + 3*v2))
            f5 = ((-3*u1 + 3*u2 - u3)*y + x*(3*v1 - 3*v2 + v3))
            f6 = (9*u2*v1 - 6*u3*v1 - 9*u1*v2 + 3*u3*v2 + 6*u1*v3 - 3*u2*v3)
            f7 = -((6*u1 - 3*u2)*y + x*(-6*v1 + 3*v2))^2
            f8 = (-3*u1*y + 3*x*v1)
            f9 = ((-3*u1 + 3*u2 - u3)*y + 9*u2*v1 - 9*u1*v2 + x*(3*v1 - 3*v2 + v3))
            f10 = ((-3*u1 + 3*u2 - u3)*y + x*(3*v1 - 3*v2 + v3)) 
            f11 = ((6*u1 - 3*u2)*y + x*(-6*v1 + 3*v2)) 
            f12 = (-9*u2*v1 + 3*u3*v1 + 9*u1*v2 - 3*u1*v3)
            f13 = ((-3*u1 + 3*u2 - u3)*y + x*(3*v1 - 3*v2 + v3)) 
            f14 = ((-3*u1 + 3*u2 - u3)*y + 9*u2*v1 - 9*u1*v2 + x*(3*v1 - 3*v2 + v3))

            local eval = f1*(f2*f3 - f4*f5) + f6*(f7 + f8*f9) + f10*(f11*f12 - f13*f14)

            eval = eval * imp_sign

            return eval
        end
    else
        -- Seems like it is possible to cut a non-degenerate cubic in a way that
        -- the segment is degenerated. Is there a better way to solve this (instead
        -- of repeating code) ?
        imp = holder[n].diagonal
    end

    holder[n].implicit = imp

    -- Store transformed control points
    holder[n].x0, holder[n].y0 = u0, v0
    holder[n].x1, holder[n].y1 = u1, v1
    holder[n].x2, holder[n].y2 = u2, v2
    holder[n].x3, holder[n].y3 = u3, v3
    holder[n].dysign = sign(v3 - v0)
end

function prepare_table.push_functions.rational_quadratic_segment(u0, v0, u1, v1, u2, v2, w, holder)
    local n = #holder + 1
    holder[n] = {}
    holder[n].type = "rational_segment"
    
    -- Compute bounding box without transforming vertices
    holder[n].xmax, holder[n].xmin = max(u0, u2), min(u0, u2)
    holder[n].ymax, holder[n].ymin = max(v0, v2), min(v0, v2)
    holder[n].dysign = sign(v2 - v0)

    -- Translate to origin
    local trans = xform.translate(-u0, -v0)
    u0, v0 = transform_point(u0, v0, trans)
    u1, v1, w = trans : apply(u1, v1, w)
    u2, v2 = transform_point(u2, v2, trans)

    holder[n].scene2origin = trans

    -- Triangle test: compute the diagonal cutting
    -- the bounding box in two triangles
    local diag_a, diag_b, diag_c = compute_implicit_line(u0, v0, u2, v2)
    holder[n].diagonal = function(x, y)
        return sign( diag_a*x + diag_b*y + diag_c )
    end

    holder[n].mid_point_diagonal = holder[n].diagonal( bezier.at2rc(0.5, u0, v0, u1, v1, w, u2, v2) )

    -- Compute implicit equation
    local det = xform.xform(u0, v0, 1, u1, v1, w, u2, v2, 1) : det()

    if det ~= 0 then
        local imp_sign = sign( 2*v2*(-u2*v1 + u1*v2) )
        holder[n].implicit = function(x, y)

            local eval = y*((4*u1^2 - 4*w*u1*u2 + u2^2)*y + 4*u1*u2*v1 - v2*4*u1^2)
            eval = eval + x*(-4*u2*v1^2 + 4*u1*v1*v2 + 
                y*(-8*u1*v1 + 4*w*u2*v1 + 4*w*u1*v2 - 2*u2*v2) + x*(4*v1^2 - 4*w*v1*v2 + v2^2))

            if imp_sign < 0 then eval = -eval end

            return eval
        end

    else
        holder[n].implicit = holder[n].diagonal
    end

    -- Store transformed control points
    holder[n].x0, holder[n].y0 = u0, v0
    holder[n].x1, holder[n].y1 = u1, v1
    holder[n].x2, holder[n].y2 = u2, v2
    holder[n].w = w
end

--------------------------------
--------- Instructions ---------
--------------------------------
function prepare_table.instructions.begin_closed_contour(shape, offset, iadd)
    local xf, data = shape.xf, shape.data
    data[offset+1], data[offset+2] = transform_point(data[offset+1], data[offset+2], xf)
end

function prepare_table.instructions.end_closed_contour(shape, offset, iadd)
    -- Fetch first vertice and then add closing edge
    local data = shape.data
    local x, y = data[offset], data[offset+1]

    local instr_offset = data[offset + 2]
    local closing_instruction = shape.offsets[iadd - instr_offset]
    local first_x, first_y = data[closing_instruction+1], data[closing_instruction+2]

    prepare_table.push_functions.linear_segment(x, y, first_x, first_y, shape.primitives, false)
end

function prepare_table.instructions.begin_open_contour(shape, offset, iadd)
    prepare_table.instructions.begin_closed_contour(shape, offset, iadd)
end

function prepare_table.instructions.end_open_contour(shape, offset, iadd)
    -- Do the same as end_closed_contour, but mark last edge as Virtual
    prepare_table.instructions.end_closed_contour(shape, offset, iadd)
    shape.primitives[ #shape.primitives ].virtual = true
end

function prepare_table.instructions.linear_segment(shape, offset, iadd)
    local data = shape.data

    -- If everything went well, first point was already transformed
    -- (by a begin_xxxx_contour instruction, for example)
    data[offset+2], data[offset+3] = transform_point(data[offset+2], data[offset+3], shape.xf)

    local x0, y0 = data[offset], data[offset+1]
    local x1, y1 = data[offset+2], data[offset+3]

    prepare_table.push_functions.linear_segment(x0, y0, x1, y1, shape.primitives, false)
end

function prepare_table.instructions.degenerate_segment(shape, offset, iadd)
    local primitives, data = shape.primitives, shape.data
    local x0, y0 = data[offset], data[offset+1]

    -- The last two parameters in data are repeated but were not transformed!!!
    data[offset+6], data[offset+7] = x0, y0

    prepare_table.push_functions.linear_segment(x0, y0, x0, y0, shape.primitives, false)
end

function prepare_table.instructions.quadratic_segment(shape, offset, iadd)
    local primitives, data = shape.primitives, shape.data

    data[offset+2], data[offset+3] = transform_point(data[offset+2], data[offset+3], shape.xf)
    data[offset+4], data[offset+5] = transform_point(data[offset+4], data[offset+5], shape.xf)

    local x0, y0 = data[offset], data[offset+1]
    local x1, y1 = data[offset+2], data[offset+3]
    local x2, y2 = data[offset+4], data[offset+5]

    -- Calculate maxima points
    local t = {}
    t[1], t[4] = 0, 1
    t[2] = truncate_parameter( (x0-x1) / (x0 - 2*x1 + x2) )
    t[3] = truncate_parameter( (y0-y1) / (y0 - 2*y1 + y2) )
    table.sort( t )

    -- Split bézier
    for i = 2, 4 do
        if t[i-1] ~= t[i] then
            u0, v0, u1, v1, u2, v2 = bezier.cut2(t[i-1], t[i], x0, y0, x1, y1, x2, y2)
            prepare_table.push_functions.quadratic_segment(u0, v0, u1, v1, u2, v2, primitives)
        end
    end
end

function prepare_table.instructions.cubic_segment(shape, offset, iadd)
    local primitives, data = shape.primitives, shape.data

    for i = 2, 6, 2 do
        data[offset+i], data[offset+(i+1)] = transform_point(data[offset+i], data[offset+(i+1)], shape.xf)        
    end

    local x0, y0 = data[offset], data[offset+1]
    local x1, y1 = data[offset+2], data[offset+3]
    local x2, y2 = data[offset+4], data[offset+5]
    local x3, y3 = data[offset+6], data[offset+7]

    -- Cross product of control points
    local d2 = 3*( x3*(2*y1 - y2 - y0) + x2*(2*y0 - 3*y1 + y3) + x0*(y1 - 2*y2 + y3) - x1*(y0 - 3*y2 + 2*y3) )
    local d3 = 3*( (3*x2 - x3)*(y0 - y1) - x1*(2*y0 - 3*y2 + y3) + x0*(2*y1 - 3*y2 + y3) )
    local d4 = 9*( x2*(y0 - y1) + x0*(y1 - y2) + x1*(y2 - y0) )

    -- Watchout for degenerated cubics! (What about degenerating into a point?
    -- just ignoring it for a while)
    if d2 == 0 and d3 == 0 then
        if d4 ~= 0 then
            -- Cubic degenerates into a quadratic (this code chunk should be
            -- merged with the one in push_quadratic_segment after)

            local u0, v0, u1, v1, u2, v2
            u0, v0, u2, v2 = x0, y0, x3, y3

            -- intersection of tangents
            u1 = (x0*(x3*(-y1 + y2) + x2*(y1 - y3)) + x1*(x3*(y0 - y2) + x2*(-y0 + y3)))
            u1 = u1 / (-(x2 - x3)*(y0 - y1) + (x0 - x1)*(y2 - y3))

            v1 = (x3*(y0 - y1)*y2 + x0*y1*y2 - x2*y0*y3 - x0*y1*y3 + x2*y1*y3 + x1*y0*(-y2 + y3))
            v1 = v1 / (-(x2 - x3)*(y0 - y1) + (x0 - x1)*(y2 - y3))

            -- Calculate maxima points
            local t = {}
            t[1], t[4] = 0, 1
            t[2] = truncate_parameter( (u0-u1) / (u0 - 2*u1 + u2) )
            t[3] = truncate_parameter( (v0-v1) / (v0 - 2*v1 + v2) )
            table.sort( t )

            -- Split bézier
            for i = 2, 4 do
                if t[i-1] ~= t[i] then
                    local a, b, c, d, e, f = bezier.cut2(t[i-1], t[i], u0, v0, u1, v1, u2, v2)
                    prepare_table.push_functions.quadratic_segment(a, b, c, d, e, f, primitives)
                end
            end

            return

        else
            -- Cubic degenerates into a line
            prepare_table.push_functions.linear_segment(x0, y0, x3, y3, shape.primitives, true)
            return
        end
    end

    -- Calculate maxima, double points and inflections
    local t = {}
    t[2], t[3] = compute_cubic_maxima(x0, x1, x2, x3)
    t[4], t[5] = compute_cubic_maxima(y0, y1, y2, y3)
    t[6], t[7] = compute_cubic_inflections(x0, y0, x1, y1, x2, y2, x3, y3, d2, d3, d4)
    t[8], t[9] = compute_cubic_doublepoint(x0, y0, x1, y1, x2, y2, x3, y3, d2, d3, d4) 
    t[1], t[10] = 0, 1

    table.sort( t )

    for i = 2, #t do
        if t[i-1] ~= t[i] then
            u0, v0, u1, v1, u2, v2, u3, v3 = bezier.cut3(t[i-1], t[i], x0, y0, x1, y1, x2, y2, x3, y3)
            prepare_table.push_functions.cubic_segment(u0, v0, u1, v1, u2, v2, u3, v3, primitives)
        end
    end
end

function prepare_table.instructions.rational_quadratic_segment(shape, start, iadd)
    local primitives, data, xf = shape.primitives, shape.data, shape.xf

    -- Unpack and transform values
    data[start+2], data[start+3] = transform_point(data[start+2], data[start+3], xf)
    data[start+5], data[start+6] = transform_point(data[start+5], data[start+6], xf)

    local x0, y0 = data[start], data[start+1]
    local x1, y1 = data[start+2], data[start+3]
    local w = data[start+4]
    local x2, y2 = data[start+5], data[start+6]

    -- Find maxima
    t = {}
    t[1], t[6] = 0, 1
    t[2], t[3] = compute_rational_maxima(x0, x1, x2, w)
    t[4], t[5] = compute_rational_maxima(y0, y1, y2, w)
    table.sort( t )

    for i = 2, 6 do
        if t[i-1] ~= t[i] then
            local u0, v0, u1, v1, r, u2, v2 = bezier.cut2rc(t[i-1], t[i], x0, y0, x1, y1, w, x2, y2)
            prepare_table.push_functions.rational_quadratic_segment(u0, v0, u1, v1, u2, v2, r, primitives)            
        end
    end
end

-----------------------------------------------------------------------------------------
-------------------------------- GEOMETRY PREPROCESSING ---------------------------------
-----------------------------------------------------------------------------------------
function prepare_table.triangle(element)
    -- We can a transformation that maps to a canonical 
    --triangle for speed!
    local shape = element.shape

    -- Transform vertices
    local x0, y0 = shape.x1, shape.y1
    local x1, y1 = shape.x2, shape.y2
    local x2, y2 = shape.x3, shape.y3

    x0, y0 = transform_point(x0, y0, shape.xf)
    x1, y1 = transform_point(x1, y1, shape.xf)
    x2, y2 = transform_point(x2, y2, shape.xf)

    -- Precompute implicit edges    
    shape.implicit = {}
    
    local compute_implicit = function(x0, y0, x1, y1)
        
        local a, b= y1-y0, -(x1-x0)
        local c = -a*x0-b*y0

        local n = #shape.implicit+1
        shape.implicit[n] = {}
        shape.implicit[n].a = a
        shape.implicit[n].b = b
        shape.implicit[n].c = c
    end

    compute_implicit(x0, y0, x1, y1)
    compute_implicit(x1, y1, x2, y2)
    compute_implicit(x2, y2, x0, y0)

    -- Bounding box info
    shape.xmax, shape.xmin = max(x2, max(x1, x0)), min(x2, min(x1, x0))
    shape.ymax, shape.ymin = max(y2, max(y1, y0)), min(y2, min(y1, y0))
end

function prepare_table.circle(element)
    -- Precompute inverse (we could precompute a transformation
    -- which maps to a canonical circle, also)
    local shape = element.shape
    shape.inversexf = shape.xf : inverse()
end

function prepare_table.path(element)
    local shape = element.shape
    shape.primitives = {}

    -- Build primitives
    for i, v in ipairs(shape.instructions) do
        local offset = shape.offsets[i]
        prepare_table.instructions[v](shape, offset, i)
    end
end

function prepare_table.polygon(element)
    local shape, data = element.shape, element.shape.data
    shape.primitives = {}

    data[1], data[2] = transform_point(data[1], data[2], shape.xf)

    -- We build primitives just as for linear_segments. We'll "trick" the
    -- rasterization function for paths.
    for i = 3, #shape.data, 2 do
        data[i], data[i+1] = transform_point(data[i], data[i+1], shape.xf)

        local x0, y0 = data[i-2], data[i-1]
        local x1, y1 = data[i], data[i+1]

        prepare_table.push_functions.linear_segment(x0, y0, x1, y1, shape.primitives, true)
    end

    -- Push closing edge
    prepare_table.push_functions.linear_segment(data[#data-1], data[#data], 
                                            data[1], data[2], shape.primitives, true)
end

-----------------------------------------------------------------------------------------
-------------------------------- PAINT PREPROCESSING ------------------------------------
-----------------------------------------------------------------------------------------
prepare_table.prepare_paint = {}

function prepare_table.prepare_paint.solid(paint, scenexf)
    -- Just multiply color alpha channel by layer opacity
    paint.data[4] = paint.data[4] * paint.opacity
end

function prepare_table.prepare_paint.lineargradient(paint, scenexf)
    local data = paint.data
    local p1, p2 = data.p1, data.p2

    fix_ramp( data.ramp )

    -- Translate p1 to the origin
    local trans = xform.translate(-p1[1], -p1[2])
    data.p1[1], data.p1[2] = transform_point(p1[1], p1[2], trans)
    data.p2[1], data.p2[2] = transform_point(p2[1], p2[2], trans)

    -- Rotate
    local rot = xform.identity()
    local dist_center = math.sqrt( p2[1]^2 + p2[2]^2 )

    if dist_center ~= 0 then
        local cos_theta = p2[1] / dist_center
        local sin_theta = math.sqrt( 1 - cos_theta^2 )

        -- Rotate clockwise if center is above <1, 0>, 
        -- counter-clockwise otherwise
        if p2[2] > 0 then sin_theta = -sin_theta end

        rot = xform.rotate( cos_theta, sin_theta )
    end

    data.p1[1], data.p1[2] = transform_point(p1[1], p1[2], rot)
    data.p2[1], data.p2[2] = transform_point(p2[1], p2[2], rot)    

    -- Precompute values
    data.grad_length = p2[1] - p1[1]

    -- "Undo" shape transformation and precompute inverse
    local canonize = rot * trans
    data.scene_to_grad = canonize * paint.xf : inverse() * scenexf : inverse()
end

function prepare_table.prepare_paint.radialgradient(paint, scenexf)
    local data, to_grad = paint.data, paint.xf
    local ramp, c, f, r = data.ramp, data.center, data.focus, data.radius

    fix_ramp( data.ramp )

    -- Translate center to the origin
    local trans = xform.translate(-f[1], -f[2])
    c[1], c[2] = transform_point(c[1], c[2], trans)
    f[1], f[2] = transform_point(f[1], f[2], trans)

    -- Compute angle between transformed focus and x axis.
    -- If focus and center coincide, do not rotate.
    local rot = xform.identity()
    local dist_center = math.sqrt( c[1]^2 + c[2]^2 )

    if dist_center ~= 0 then
        local cos_theta = c[1] / dist_center
        local sin_theta = math.sqrt( 1 - cos_theta^2 )

        -- Rotate clockwise if center is above <1, 0>, 
        -- counter-clockwise otherwise
        if c[2] > 0 then sin_theta = -sin_theta end

        rot = xform.rotate( cos_theta, sin_theta )
    end

    c[1], c[2] = transform_point(c[1], c[2], rot)
    f[1], f[2] = transform_point(f[1], f[2], rot)

    -- Store transform and its inverse. We'll transform the point using the
    -- "direct" one, and we'll use the inverse to compose with other transformations
    local canonize = rot * trans
    data.scene_to_grad = canonize * (to_grad : inverse()) * (scenexf : inverse())
end

function prepare_table.prepare_paint.texture(paint, scenexf)
    local data, tex2scene = paint.data, paint.xf

    paint.scene2tex = (tex2scene : inverse() * scenexf : inverse())
end

-- prepare scene for sampling and return modified scene
local function preparescene(scene)

    -- Erase this after
    export_cell(0)
    -- Erase this after

    for i, element in ipairs(scene.elements) do
        element.shape.xf = scene.xf * element.shape.xf
        prepare_table[element.shape.type](element)
        prepare_table.prepare_paint[element.paint.type](element.paint, scene.xf)
    end

    return scene
end

-----------------------------------------------------------------------------------------
-------------------------------- PATH SAMPLING ------------------------------------------
-----------------------------------------------------------------------------------------
local sample_table = {}
sample_table.sample_path = {}

-- TODO: MUST FUSION ALL THESE TESTS AFTER

function sample_table.sample_path.linear_segment(primitive, x, y)
    -- Bounding box tests
    if y >= primitive.ymax or y < primitive.ymin then return 0 end
    if x > primitive.xmax then return 0 end
    if x <= primitive.xmin then return primitive.dysign end

    -- Implicit test
    local eval = sign( primitive.a * x + primitive.b * y + primitive.c )

    if eval < 0 then return primitive.dysign
    else return 0 end
end

function sample_table.sample_path.degenerate_segment(primitive, x, y)
    if y == primitive.y and x <= primitive.x then
        return sign (primitive.dy0 )
    else
        return 0
    end
end

-- TODO: Fusion these two functions
function sample_table.sample_path.quadratic_segment(primitive, x, y)

    -- Bounding box test
    if y >= primitive.ymax or y < primitive.ymin then return 0 end
    if x > primitive.xmax then return 0 end
    if x <= primitive.xmin then return primitive.dysign end

    local x0, y0 = primitive.x0, primitive.y0
    local x1, y1 = primitive.x1, primitive.y1
    local x2, y2 = primitive.x2, primitive.y2

    x, y = transform_point(x, y, primitive.scene_to_canonic)

    -- Triangle test -> skip if point is inside the triangle fully covered
    -- (or fully uncovered) by Bézier
    local point_diagonal = primitive.diagonal(x, y)
    if point_diagonal == -primitive.mid_point_diagonal then
        if point_diagonal < 0 then return primitive.dysign
        else return 0 end
    end

    -- Implicit test
    local eval = primitive.implicit(x, y)

    if eval < 0 then
        return primitive.dysign
    else return 0 end
end

function sample_table.sample_path.cubic_segment(primitive, x, y)
    if y >= primitive.ymax or y < primitive.ymin then return 0 end
    if x > primitive.xmax then return 0 end
    if x <= primitive.xmin then return primitive.dysign end

    x, y = transform_point(x, y, primitive.scene_to_canonic)

    local x0, y0 = primitive.x0, primitive.y0
    local x1, y1 = primitive.x1, primitive.y1
    local x2, y2 = primitive.x2, primitive.y2
    local x3, y3 = primitive.x3, primitive.y3

    -- First triangle test: skip if point is inside the triangle fully
    -- covered (or fully uncovered) by Bézier
    local point_diagonal = primitive.diagonal(x, y)
    
    if point_diagonal == -primitive.mid_point_diagonal then
        if point_diagonal < 0 then return primitive.dysign
        else return 0 end
    end

    -- Second triangle test: if point is outside it, return 0 or dysign
    -- depending whether point is to the left or to the reight to the curve;
    -- otherwise, evaluate implicitly
    if primitive.inside_triangle(x, y) == false then
        if point_diagonal < 0 then 
            return primitive.dysign
        elseif point_diagonal > 0 then
            return 0 
        else
            -- Point is exactly on the diagonal. In this case, we return 0 or +-1
            -- depending on the position of point t = 0.5
            if primitive.mid_point_diagonal > 0 then return primitive.dysign
            else return 0 end
        end
    else
        local eval = primitive.implicit(x, y)
        if eval > 0 then return primitive.dysign
        end return 0
    end
end

function sample_table.sample_path.rational_segment(primitive, x, y)

    -- This function is essentialy the same as quadratic_segment(), so we
    -- could fusion both
    if y >= primitive.ymax or y < primitive.ymin then return 0 end
    if x > primitive.xmax then return 0 end
    if x <= primitive.xmin then return primitive.dysign end

    x, y = transform_point(x, y, primitive.scene2origin)

    -- Triangle test
    local point_diagonal = primitive.diagonal(x, y)
    if point_diagonal == -primitive.mid_point_diagonal then
        if point_diagonal < 0 then return primitive.dysign
        else return 0 end
    end

    -- Implicit test
    if primitive.implicit(x, y) < 0 then
        return primitive.dysign
    else return 0 end
end

-----------------------------------------------------------------------------------------
------------------------------- PAINT SAMPLING ------------------------------------------
-----------------------------------------------------------------------------------------
sample_table.sample_paint = {}
sample_table.sample_paint.spread_table = {}

sample_table.sample_paint.spread_table = {
    ["repeat"] = function(v)
        if v > 1 then return v - floor(v) 
        elseif v < 0 then return 1 + (v - ceil(v))
        else return v end
    end
}

function sample_table.sample_paint.spread_table.reflect(v)
    if v >= 0 then
        local int = floor(v)
        if int % 2 == 0 then return v - int
        else return 1 - (v - int) end
    else
        local int = ceil(v)
        if int % 2 == 0 then return -(v - int)
        else return 1-(v - int) end
    end
end

function sample_table.sample_paint.spread_table.pad(v)
    if v > 1 then return 1
    elseif v < 0 then return 0
    else return v end
end

---------------
---------------

function sample_table.sample_paint.solid(paint, x, y)
    return paint.data
end

function sample_table.sample_paint.lineargradient(paint, x, y)

    local data, ramp = paint.data, paint.data.ramp
    local p1 = data.p1
    x, y = transform_point(x, y, data.scene_to_grad)

    -- Compute ratio
    local k = (x-p1[1]) / data.grad_length

    -- Wrap, sample, interpolate
    local wrapped = sample_table.sample_paint.spread_table[ramp.spread](k)
    local off = search_in_ramp(ramp, wrapped)

    local inter_factor =  (wrapped - ramp[off])/(ramp[off+2] - ramp[off])
    local out = interpolate_colors(ramp[off+1], ramp[off+3],  inter_factor)

    out[4] = out[4] * paint.opacity

    return out
end

function sample_table.sample_paint.radialgradient(paint, x, y)
    local data = paint.data
    local ramp, center, f, r = data.ramp, data.center, data.focus, data.radius

    x, y = transform_point(x, y, data.scene_to_grad)

    -- Compute intersection of the line passing through origin
    -- and (x,y) with the circle
    local a = x^2 + y^2
    local b = -2*(x*center[1] + y*center[2])
    local c = center[1]^2 + center[2]^2 - r^2
    local n, r1, s1, r2, s2 = quadratic.quadratic(a, b, c)

    -- We're interest in the positive root for t (the one which goes
    -- towards the point (x,y) )
    local t1, t2 = 0, 0
    if n > 0 then t1 = r1/s1 end
    if n > 1 then t2 = r2/s2 end
    local t = max(t1, t2)

    -- Ratio of [point - focus] and [inter - focus] is 1/t
    local k = 1/t

    local wrapped = sample_table.sample_paint.spread_table[ramp.spread](k)
    local off = search_in_ramp(ramp, wrapped)

    local inter_factor = (wrapped - ramp[off])/(ramp[off+2] - ramp[off]) 
    local out = interpolate_colors(ramp[off+1], ramp[off+3],  inter_factor)

    -- Compose with opacity
    out[4] = out[4] * paint.opacity

    return out
end

function sample_table.sample_paint.texture(paint, x, y)
    local data = paint.data
    local img = data.image

    x, y = transform_point(x, y, paint.scene2tex)

    local wrapped_x = sample_table.sample_paint.spread_table[data.spread](x)
    local wrapped_y = sample_table.sample_paint.spread_table[data.spread](y)

    local tex_x = wrapped_x * img.width 
    local tex_y = wrapped_y * img.height

    -- Use nearest neighbour
    local tex_x, tex_y = ceil(tex_x), ceil(tex_y)
    local r, g, b, a = img : get(tex_x, tex_y)

    a = a * paint.opacity

    return {r, g, b, a}
end

-----------------------------------------------------------------------------------------
---------------------------- GEOMETRY SAMPLING ------------------------------------------
-----------------------------------------------------------------------------------------
function sample_table.triangle(element, x, y)
    local implicit = element.shape.implicit
    local xmin, xmax = element.shape.xmin, element.shape.xmax
    local ymin, ymax = element.shape.ymin, element.shape.ymax

    -- Bounding box tests (closed bottom, open top, closed right, open left)
    if y < ymin or y >= ymax then return BGColor end
    if x <= xmin or x > xmax then return BGColor end

    -- Implicit test
    local edge_sign = {}
    for i = 1, 3 do
        edge_sign[i] = sign( implicit[i].a*x + implicit[i].b*y + implicit[i].c )
    end

    if edge_sign[1] == edge_sign[2] and edge_sign[2] == edge_sign[3] then
        return element.paint.data
    else
        return BGColor
    end
end

function sample_table.circle(element, x, y)
    local shape, paint = element.shape, element.paint
    local cx, cy, r = shape.cx, shape.cy, shape.r

    -- Map point to untransformed circle
    tx, ty = transform_point(x, y, shape.inversexf)
    local d = math.sqrt( (cx-tx)^2 + (cy-ty)^2 )

    if d <= r then
        return sample_table.sample_paint[paint.type](paint, x, y)
    else return BGColor end
end

function sample_table.path(element, x, y)
    local paint, shape, primitives = element.paint, element.shape, element.shape.primitives

    local count = 0
    for i,prim in ipairs(primitives) do
        count = count + sample_table.sample_path[prim.type](prim, x, y)
    end

    local paint_flag
    if element.type == "fill" then
        paint_flag = (count ~= 0)
    elseif element.type == "eofill" then
        paint_flag = (count % 2 ~= 0)
    end

    if paint_flag == true then
        return sample_table.sample_paint[paint.type](paint, x, y)
    else return BGColor end
end

function sample_table.polygon(element, x, y)
    return sample_table.path(element, x, y)
end

-- sample scene at x,y and return r,g,b,a
local function sample_point(scene, x, y)

    -- Deep copy of BGColor
    -- this is extremely slow and unnecessary! Must think in a way
    -- of removing this
    local out = {}
    for i, v in ipairs(BGColor) do
        out[i] = v
    end

    for i = 1, #scene.elements do
        local element = scene.elements[i]
        local temp = sample_table[element.shape.type](element, x, y)

        -- Superpose images
        if temp ~= BGColor then

            -- Alpha blend current color with new layer
            -- This uses no premultiplication! Slow
            for j = 1, 3 do 
                out[j] = alpha_composite(temp[j], temp[4], out[j], out[4])
            end

            out[4] = temp[4] + out[4]*(1 - temp[4])
        end
    end

    return out
end

-- Here we'll supersample each pixel
local function sample(scene, x, y)

    local samples = {}
    local pattern = noise[1]

    for i = 1, #pattern, 2 do
        local ind = (i+1)/2
        local dx, dy = pattern[i], pattern[i+1]

        samples[ind] = {}
        samples[ind] = sample_point(scene, x + dx, y + dy)
    end

    -- Get arithmetic average of samples
    -- do different averages/filters produce better results?
    -- e.g. gaussian filter
    local sum = {0,0,0,0}
    for i = 1, #samples do
        for j = 1, 4 do
            samples[i] = ungamma(samples[i])
            sum[j] = sum[j] + samples[i][j]
        end
    end

    for j = 1, 4 do
        sum[j] = sum[j] / #samples
    end

    sum = gamma(sum)

    return unpack(sum)
end

-----------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------
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
    -- transform and prepare scene for rendering
    scene = preparescene(scene)
    -- get viewport
    local vxmin, vymin, vxmax, vymax = unpack(viewport, 1, 4)
stderr("preprocess in %.3fs\n", time:elapsed())
time:reset()
    -- get image width and height from viewport
    local width, height = vxmax-vxmin, vymax-vymin
    -- allocate output image
    local img = image.image(width, height)
    -- render
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
