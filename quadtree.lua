local _M = {}
-- list of functions
local command = require"command"
local element = require"element"
local paint = require"paint"
local path = require"path"
local quad = require("quadratic")
local cubic = require("cubic")
local bernstein = require"bernstein"
local floor = math.floor
local lerp = bernstein.lerp1
local lerp2 = bernstein.lerp2
local lerp3 = bernstein.lerp3

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

function quadtree_meta.__index.printree(self)
    if self.children then
        for i,child in pairs(self.children) do
			local e = child:printree()
			if e~= nil then
			for i= #e,1,-1 do
						table.insert(self.elements,e[i])
			end
		end
    end
	return self.elements
    else
		local path = path.path{command.M,self.left,self.top,command.L,self.left,self.top+self.height,command.L,self.left+self.width,self.top+self.height,command.L,self.left+self.width,self.top,command.Z,command.M,self.left+1,self.top+1,command.L,self.left+1,self.top+self.height,command.L,self.left+self.width,self.top+self.height,command.L,self.left+self.width,self.top+1,command.Z}  
		local e = element.eofill(path, paint.solid({0,0,0,1}, 1))
		e.shape.xmin = self.left
		e.shape.ymin = self.top
		e.shape.xmax = self.left+self.width
		e.shape.ymax = self.top+self.height
		e.shape.parameter = {{0,0},{self.height-1,0,-(self.height-1)*(self.left+1)},{0,self.width-1,-(self.width-1)*(self.top+self.height)},{-self.height+1,0,(self.height-1)*(self.left+self.width)},{0,-self.width+1,(self.width-1)*(self.top+1)},{0,0},{self.height+1,0,-(self.height+1)*(self.left-1)},{0,self.width+1,-(self.width+1)*(self.top+self.height)},{-self.height-1,0,(self.height+1)*(self.left+self.width)},{0,-self.width-1,(self.width+1)*(self.top-1)}}
		--dump(e)
		table.insert(self.elements,e)
		--dump(scene.elements)
		return self.elements
    end
    --return scene
end

-- This is to get the y coordinates
-- works for horizontal lines
local cy = {}
cy.get  = function(a,p)
    return a[2]
end

local bt = {}
bt.get = function(u,v)
    return u > v
end

local bte = {}
bte.get = function(u,v)
    return u >= v
end

local lt = {}
lt.get = function(u,v)
    return u < v
end

local lte = {}
lte.get = function(u,v)
    return u <= v
end

local function linear_intersection(x0,x1,w)
    return (w - x0)/(x1 -x0)
end

local function quadratic_intersection(x0,x1,x2,w)
    local a,b,c = quadratic_coefficients(x0,x1,x2)
    local t
    roots = { quad.quadratic(a,b,c-w)}
    if roots[1] == 2 then
        for i = 2,4,2 do
            t = roots[i]/roots[i+1]
            if t > 0 and t < 1 then
                return t
            end
        end
    end
end

local function cubic_intersection(x0,x1,x2,x3,w)
    local a,b,c,d = cubic_coefficients(x0,x1,x2,x3)
    local t 
    roots = { cubic.cubic(a,b,c,d-w) }
    for i = 2,#roots-1,2 do 
        t = roots[i]/roots[i+1]
        if t > 0 and t < 1 then
            return t
        end
    end
    
end

local function rational_quadratic_intersection(x0,x1,w1,x2,x)
    local a,b,c = quadratic_coefficients(x0,x1,x2)
    local d,e,f = quadratic_coefficients(1,w1,1) 
    local  ca, cb, cc = a - x*d, b - x*e, c - x*f

    local root = {quad.quadratic(ca,cb,cc) }
    if root[1] == 2 then
        for i = 2,4,2 do
            local t = root[i]/root[i+1]
            if t >= 0 and t <= 1 then
                return t
            end
        end
    end
end

-- c(Cooridinate): this will be cx o cy
-- o: lt,lte,bt,bte
-- value = {xvalue,yvalue}
local function clip(c,o,value,forward)
    local fx,fy = nil,nil -- the first point inside the path
    local px,py -- the last point added to the new path

    local iterator = {}

    function iterator:begin_closed_contour(len,x0,y0)
        if o.get(c.get({x0,y0}),c.get(value)) then
            px, py = x0,y0
            fx, fy = x0,y0
            forward:begin_closed_contour(_,x0,y0)
        end
    end
    iterator.begin_open_contour = iterator.begin_closed_contour
    function iterator:linear_segment(x0,y0,x1,y1)
        if o.get(c.get({x0,y0}),c.get(value))
            and o.get(c.get({x1,y1}),c.get(value)) then
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end
            forward:linear_segment(x0,y0,x1,y1)
            px,py = x1,y1
        elseif not o.get(c.get({x0,y0}),c.get(value))
            and o.get(c.get({x1,y1}),c.get(value)) then
            local t = linear_intersection(c.get({x0,y0}),c.get({x1,y1}),c.get(value))
            local px0 = lerp(x0,x1,t)
            local py0 = lerp(y0,y1,t)
            if fx == nil then
                fx, fy = px0,py0
                forward:begin_closed_contour(_,fx,fy)
            else
            forward:linear_segment(px,py,px0,py0)
            end
            forward:linear_segment(px0,py0,x1,y1)
            px,py = x1,y1
        elseif o.get(c.get({x0,y0}),c.get(value)) 
            and not o.get(c.get({x1,y1}),c.get(value)) then
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end
            local t = linear_intersection(c.get({x0,y0}),c.get({x1,y1}),c.get(value))
            local px1 = lerp(x0,x1,t)
            local py1 = lerp(y0,y1,t)
            forward:linear_segment(x0,y0,px1,py1)
            px,py = px1,py1
        end
    end

    function iterator:quadratic_segment(x0,y0,x1,y1,x2,y2)
        if o.get(c.get({x0,y0}),c.get(value)) 
            and o.get(c.get({x2,y2}),c.get(value)) then
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end
            forward:quadratic_segment(x0,y0,x1,y1,x2,y2)
            px,py = x2,y2
        elseif not o.get(c.get({x0,y0}),c.get(value)) 
            and o.get(c.get({x2,y2}),c.get(value)) then
            local t = quadratic_intersection(c.get({x0,y0}),
                c.get({x1,y1}),c.get({x2,y2}), c.get(value))
            local px0 = lerp2(x0,x1,x2,t,t)
            local py0 = lerp2(y0,y1,y2,t,t)

            local px1 = lerp2(x0,x1,x2,t,1)
            local py1 = lerp2(y0,y1,y2,t,1)
            if fx == nil then
                fx, fy = px0,py0
                forward:begin_closed_contour(_,fx,fy)
            else
                forward:linear_segment(px,py,px0,py0)
            end
            forward:quadratic_segment(px0,py0,px1,py1,x2,y2)
            px,py = x2,y2

        elseif o.get(c.get({x0,y0}),c.get(value)) 
            and not o.get(c.get({x2,y2}),c.get(value)) then
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end
            local t = quadratic_intersection(c.get({x0,y0}),
                c.get({x1,y1}),c.get({x2,y2}), c.get(value))
            local px1 = lerp2(x0,x1,x2,0,t)
            local py1 = lerp2(y0,y1,y2,0,t)

            local px2 = lerp2(x0,x1,x2,t,t)
            local py2 = lerp2(y0,y1,y2,t,t)
            forward:quadratic_segment(x0,y0,px1,py1,px2,py2)
            px,py = px2,py2
        end
    end

    function iterator:rational_quadratic_segment(x0,y0,x1,y1,w1,x2,y2)
        if o.get(c.get({x0,y0}),c.get(value)) and o.get(c.get({x2,y2}),c.get(value)) then
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end
            forward:rational_quadratic_segment(x0,y0,x1,y1,w1,x2,y2)
            px,py = x2,y2
        elseif not o.get(c.get({x0,y0}),c.get(value)) and o.get(c.get({x2,y2}),c.get(value)) then
            local t = rational_quadratic_intersection(c.get({x0,y0}),c.get({x1,y1}),w1,c.get({x2,y2}),c.get(value)) 
            local px0,py0,px1,py1,pw1,px2,py2 = cutr2s(t,1,x0,y0,x1,y1,w1,x2,y2)
            if fx == nil then
                fx, fy = px0,py0
                forward:begin_closed_contour(_,fx,fy)
            else
                forward:linear_segment(px,py,px0,py0)
            end
            forward:rational_quadratic_segment(px0,py0,px1,py1,pw1,x2,y2)
            px,py = x2,y2
        elseif o.get(c.get({x0,y0}),c.get(value))  and not o.get(c.get({x2,y2}),c.get(value)) then
            local t = rational_quadratic_intersection(c.get({x0,y0}),c.get({x1,y1}),w1,c.get({x2,y2}),c.get(value)) 
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end
            local px0,py0,px1,py1,pw1,px2,py2 = cutr2s(0,t,x0,y0,x1,y1,w1,x2,y2)
            forward:rational_quadratic_segment(x0,y0,px1,py1,pw1,px2,py2)
            px,py = px2,py2

        end
    end

    function iterator:cubic_segment(x0,y0,x1,y1,x2,y2,x3,y3)
        if o.get(c.get({x0,y0}),c.get(value)) 
            and o.get(c.get({x3,y3}),c.get(value)) then
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end
            forward:cubic_segment(x0,y0,x1,y1,x2,y2,x3,y3)
            px,py = x3,y3
        elseif not o.get(c.get({x0,y0}),c.get(value)) 
            and o.get(c.get({x3,y3}),c.get(value)) then

            local t = cubic_intersection(c.get({x0,y0}),c.get({x1,y1}),c.get({x2,y2}),c.get({x3,y3}),c.get(value))

            local px0 = lerp3(x0,x1,x2,x3,t,t,t)
            local py0 = lerp3(y0,y1,y2,y3,t,t,t)
            if fx == nil then
                fx, fy = px0,py0
                forward:begin_closed_contour(_,fx,fy)
            else
                forward:linear_segment(px,py,px0,py0)
            end

            local px1 = lerp3(x0,x1,x2,x3,t,t,1)
            local py1 = lerp3(y0,y1,y2,y3,t,t,1)

            local px2 = lerp3(x0,x1,x2,x3,t,1,1)
            local py2 = lerp3(y0,y1,y2,y3,t,1,1)

            forward:cubic_segment(px0,py0,px1,py1,px2,py2,x3,y3)
            px,py = x3,y3
        elseif o.get(c.get({x0,y0}),c.get(value)) 
            and not o.get(c.get({x3,y3}),c.get(value)) then
            if fx == nil then
                fx, fy = x0,y0
                forward:begin_closed_contour(_,fx,fy)
            end

            local t = cubic_intersection(c.get({x0,y0}),c.get({x1,y1}),c.get({x2,y2}),c.get({x3,y3}),c.get(value))

            local px1 = lerp3(x0,x1,x2,x3,0,0,t)
            local py1 = lerp3(y0,y1,y2,y3,0,0,t)

            local px2 = lerp3(x0,x1,x2,x3,0,t,t)
            local py2 = lerp3(y0,y1,y2,y3,0,t,t)

            local px3 = lerp3(x0,x1,x2,x3,t,t,t)
            local py3 = lerp3(y0,y1,y2,y3,t,t,t)

            forward:cubic_segment(x0,x0,px1,py1,px2,py2,px3,py3)
            px,py = px3,py3
        end
    end

    function iterator:end_closed_contour(len)
        if px ~= fx or py ~= fy then
            forward:linear_segment(px,py,fx,fy)
        end
        forward:end_closed_contour(_)
    end
    iterator.end_open_contour = iterator.end_closed_contour
    return iterator
end

local function clippath(c,o,value,oldpath)
    local newpath = _M.path()
    newpath:open()
    oldpath:iterate(
        clip(c,o,value,
            newcleaner(newpath)
            )
    )
    newpath:close()
    return newpath
end

function quadtree_meta.__index.scenetoleaf(self,scene, vxmin, vymin, vwidth, vheight,maxdepth)
	local depth = 0	
	if maxdepth >= 0 then
		for k = #scene.elements,1,-1 do
			local e =  scene.elements[k]  -- Verifica se tem elemento dentro da Folha
    			if  e.shape.xmax > vxmin and vxmin+vwidth > e.shape.xmin  and  e.shape.ymax > vymin and vymin+vheight > e.shape.ymin then
				self:subdivide(depth,maxdepth)
				--dump(e) 
				table.insert(self.elements,e)
			end
		end
		-- Verifica se tem 
		if self.children ~= nil then
			for k = #scene.elements,1,-1 do
			local e =  self.elements[k]
			for i,child in pairs(self.children) do
				local path
				print(child.left,child.top)
				if i == 1 then
					path = clippath(cx,bt,{child.left,child.top},e.shape)
				--	table.insert(quadtree.elements,e)
				elseif i == 2 then
					path = clippath(cx,bt,{child.left,child.top},e.shape)
				--	table.insert(quadtree.elements,e)
				elseif i == 3 then
					path = clippath(cx,bt,{child.left,child.top},e.shape)
				--	table.insert(quadtree.elements,e)
				else 
					path = clippath(cx,bt,{child.left,child.top},e.shape)
				--	table.insert(quadtree.elements,e)
				end
				--local path = clippath(cx,lt,{child.left,child.top},e.shape)
				dump(path)				
				path.xmin = max(e.shape.xmin,child.left) 
				path.xmax = min(e.shape.xmax,child.left+child.width)
				path.ymin = max(e.shape.xmin,child.top)
				path.ymax = min(e.shape.xmax,child.top+child.height)
				dump(path)
				path.parameter = e.shape.parameter
				table.insert(self.children[i].elements,element.fill(path,e.paint))
				--print("Passei Aqui")
				self.children[i] = scenetoleaf(quadtree,child.left,child.top,child.width,child.height, maxdepth-1)
			end
			end
		end
	end
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
