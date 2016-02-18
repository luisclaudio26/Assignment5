local path = require"path"

local _M = {}

--------------------------------------------------------------------------
--------------------- "CONVERT" PRIMITIVE INTO PATH ----------------------
--------------------------------------------------------------------------
local primitive2path = {}
function primitive2path.linear_segment(prim)
	local x0, y0, x1, y1 = prim.x0, prim.y0, prim.x1, prim.y1

	local built = path.path()

	built : begin_open_contour(0, x0, y0)
	built : linear_segment(x0,y0,x1,y1)
	built : end_open_contour(x1,y1,0)

	built = built : stroke(4)

	return built
end

function primitive2path.quadratic_segment(prim)
    local built = path.path()

    built : begin_open_contour(0, prim.x0, prim.y0)
    built : quadratic_segment(prim.x0, prim.y0, prim.x1, prim.y1, prim.x2, prim.y2)
    built : end_open_contour(prim.x2, prim.y2, 0)

    built = built : stroke(3)

    return built
end

function primitive2path.cubic_segment(prim)
   local built = path.path()

    built : begin_open_contour(0, prim.x0, prim.y0)
    built : cubic_segment(prim.x0, prim.y0, prim.x1, prim.y1, prim.x2, prim.y2, prim.x3, prim.y3)
    built : end_open_contour(prim.x3, prim.y3, 0)

    built = built : stroke(3)

    return built 
end


function primitive2path.rational_quadratic_segment(prim)
    local built = path.path()

    built : begin_open_contour(0, prim.x0, prim.y0)
    built : rational_quadratic_segment(prim.x0, prim.y0, prim.x1, prim.y1, prim.w, prim.x2, prim.y2)
    built : end_open_contour(prim.x2, prim.y2, 0)

    built = built : stroke(3)

    return built
end

------------------------------
------------------------------
function _M.export_cell(cell)
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

    local complete_scene = scene.scene(elements)
    local output = assert(io.open("cell.svg", "wb"))

    require("svg").render(complete_scene, {cell.xmin, cell.ymin, cell.xmax, cell.ymax}, output)

    output.close()
end

return _M