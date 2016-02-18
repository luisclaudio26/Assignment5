local path = require"path"

local _M = {}


local primitive2path = {}
function primitive2path.linear_segment(prim)
	local x0, y0, x1, y1 = prim.x0, prim.y0, prim.x1, prim.y1

	local built = path.path()

	built : begin_open_contour(0, x0, y0)
	built : linear_segment(x0,y0,x1,y1)
	built : end_open_contour(x1,y1,0)

	built = built : stroke(4)

	--print(#built.instructions)
	return built

end



function _M.export_cell(cell)

	local element, scene = require("element"), require("scene")

	local elements = {}
    for i, v in ipairs(cell.shapes) do

        -- Quero um ELEMENT com o SEGMENT e pintura tipo STROKE
        -- Antes disso, quero gerar um PATH a partir de uma PRIMITIVE qualquer

        -- ISTO ASSUME QUE v.segment É UMA PRIMITIVE DA MANEIRA QUE CONSTRUO
        -- NAS FUNÇÕES DA TABELA PUSH_FUNCTIONS
        local p = primitive2path[v.segment.type]( v.segment )

        if v.fill_type == "fill" then
	        elements[i] = element.fill(p, v.paint)
    	end
    end

    local complete_scene = scene.scene(elements)

    local output = assert(io.open("cell.rvg", "wb"))

    require("rvg").render(complete_scene, {cell.xmin, cell.ymin, cell.xmax, cell.ymax}, output)

    output.close()
end

return _M