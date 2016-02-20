local _M = {}
-- list of functions
local command = require"command"
local element = require"element"
local paint = require"paint"
local path = require"path"
local floor = math.floor

local quadtree_meta = {}
quadtree_meta.__index = {}


function _M.quadtree(_left, _top, _width, _height)
    return setmetatable({
            left = _left,
            top = _top,
            width = _width,
            height = _height,
            children = nil,
            objects = {},
            elements = {}
    }, quadtree_meta)
end




function quadtree_meta.__index.subdivide(self,depth,maxdepth)
    if self.children == nil and depth < maxdepth then
        local x = self.left
        local y = self.top
        local w = floor(self.width / 2)
        local h = floor(self.height / 2)
        self.children = {
            _M.quadtree(x , y , w, h),
            _M.quadtree(x + w, y , w, h),
            _M.quadtree(x , y + h, w, h),
            _M.quadtree(x + w, y + h, w, h)
        }
    end
end
-- verifica se tem o ponto e retonar o quadtree que o possui.
function quadtree_meta.__index.check(self,x,y)
    if self.children then
		for i,child in pairs(self.children) do
			local chi = child:check(x,y)
			if chi ~= nil then 
				return chi
			end
		end
    else
		if x-self.left < self.width and y - self.top  < self.height then 
			return self
		else
			return nil
		end
    end
end

function quadtree_meta.__index.printree(self,scene)
    if self.children then
        for i,child in pairs(self.children) do
	    local e = child:printree(self,scene)
	    if e~= nil then
		for i= #e,1,-1 do
            		table.insert(scene.elements,e[i])
		end
	    end
	    --dump(scene.elements		
        end
	return scene.elements
    else
	local path = path.path{command.M,self.left,self.top,command.L,self.left,self.top+self.height,command.L,self.left+self.width,self.top+self.height,command.L,self.left+self.width,self.top,command.Z,command.M,self.left+1,self.top+1,command.L,self.left+1,self.top+self.height,command.L,self.left+self.width,self.top+self.height,command.L,self.left+self.width,self.top+1,command.Z}  
	local e = element.eofill(path, paint.solid({0,0,0,1}, 1))
	e.shape.xmin = self.left
	e.shape.ymin = self.top
	e.shape.xmax = self.left+self.width
	e.shape.ymax = self.top+self.height
	e.shape.parameter = {{0,0},{self.height-1,0,-(self.height-1)*(self.left+1)},{0,self.width-1,-(self.width-1)*(self.top+self.height)},{-self.height+1,0,(self.height-1)*(self.left+self.width)},{0,-self.width+1,(self.width-1)*(self.top+1)},{0,0},{self.height+1,0,-(self.height+1)*(self.left-1)},{0,self.width+1,-(self.width+1)*(self.top+self.height)},{-self.height-1,0,(self.height+1)*(self.left+self.width)},{0,-self.width-1,(self.width+1)*(self.top-1)}}
	--dump(e)
	table.insert(scene.elements,e)
	--dump(scene.elements)
	return scene.elements
    end
    --return scene
end

function quadtree_meta.__index.clean() 
    if self.children and depth < maxdepth then
        for i,child in pairs(self.children) do
            child:clean()
        end
    else
        if #self.element == 0 then
            self = nil
        end
    end
end


return _M
