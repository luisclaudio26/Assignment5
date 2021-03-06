local path = require"path"

local _M = {}

--------------------------------------------------------------------------
--------------------- "CONVERT" PRIMITIVE INTO PATH ----------------------
--------------------------------------------------------------------------
local primitive2path = {}
function primitive2path.linear_segment(prim, thick)
    thick = thick or 3

	-- local x0, y0, x1, y1 = prim.x0, prim.y0, prim.x1, prim.y1
    -- Adaptação pro código da Alana:
    local x0, x1 = table.unpack(prim.x)
    local y0, y1 = table.unpack(prim.y)

	local built = path.path()

	built : begin_open_contour(0, x0, y0)
	built : linear_segment(x0,y0,x1,y1)
	built : end_open_contour(x1,y1,0)

	built = built : stroke(thick)

	return built
end

function primitive2path.quadratic_segment(prim, thick)
    thick = thick or 3
    local built = path.path()

    built : begin_open_contour(0, prim.x0, prim.y0)
    built : quadratic_segment(prim.x0, prim.y0, prim.x1, prim.y1, prim.x2, prim.y2)
    built : end_open_contour(prim.x2, prim.y2, 0)

    built = built : stroke(thick)

    return built
end

function primitive2path.cubic_segment(prim, thick)
    thick = thick or 3
    local built = path.path()

    built : begin_open_contour(0, prim.x0, prim.y0)
    built : cubic_segment(prim.x0, prim.y0, prim.x1, prim.y1, prim.x2, prim.y2, prim.x3, prim.y3)
    built : end_open_contour(prim.x3, prim.y3, 0)

    built = built : stroke(thick)

    return built 
end


function primitive2path.rational_quadratic_segment(prim, thick)
    thick = thick or 3
    local built = path.path()

    built : begin_open_contour(0, prim.x0, prim.y0)
    built : rational_quadratic_segment(prim.x0, prim.y0, prim.x1, prim.y1, prim.w, prim.x2, prim.y2)
    built : end_open_contour(prim.x2, prim.y2, 0)

    built = built : stroke(thick)

    return built
end

------------------------------
------------------------------
local function drawBoundingBox(xmin, ymin, xmax, ymax, elements)
    local paint, color, element = require("paint"), require("color"), require("element")

    -- Bounding box
    local l = {}
    l[1] = primitive2path.linear_segment({x0 = xmin, y0 = ymin, x1 = xmax, y1 = ymin},0.7)
    l[2] = primitive2path.linear_segment({x0 = xmax, y0 = ymin, x1 = xmax, y1 = ymax},0.7)
    l[3] = primitive2path.linear_segment({x0 = xmax, y0 = ymax, x1 = xmin, y1 = ymax},0.7)
    l[4] = primitive2path.linear_segment({x0 = xmin, y0 = ymax, x1 = xmin, y1 = ymin},0.7)

    for i,v in ipairs(l) do
        elements[#elements + 1] = element.fill(v, paint.solid(color.rgb8(255,0,0)) )
    end
end

function _M.export_cell(cell, filename)

    filename = filename or "cell"

	local element, scene = require("element"), require("scene")

	local elements = {}
    for i, v in ipairs(cell.shapes) do
        -- [Assumption] Assumes v.segment is a primitive built by some function
        -- inside prepare_table.push_functions. Obviously, this can be adapted
        -- by changing the functions inside primitive2path to read the correct
        -- values.
        local p = primitive2path[v.segment.type]( v.segment )

        if v.fill_type == "fill" then
	        elements[i] = element.fill(p, v.paint)
        elseif v.fill_type == "eofill" then
            elements[i] = element.eofill(p, v.paint)
    	end
    end

    drawBoundingBox(cell.xmin, cell.ymin, cell.xmax, cell.ymax, elements)

    local complete_scene = scene.scene(elements)
    local output = assert(io.open("cells\\" .. filename .. ".svg", "wb"))

    require("svg").render(complete_scene, {cell.xmin, cell.ymin, cell.xmax, cell.ymax}, output)

    output.close()
end

return _M
