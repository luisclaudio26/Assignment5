local driver = require"driver"
local image = require"image"
local chronos = require"chronos"
local xform = require("xform")
local vector = require"vector"
local bezier = require"bezier"
local quad = require("quadratic")
local cubic = require("cubic")
local bernstein = require"bernstein"
local util = require"util"
local blue = require"blue"
local filter = require"filter"

local samples =  blue[1] -- Mude o sampling alterando de blue[{1,8,16,32}]
local power3 = bernstein.power3
local stderr = util.stderr
local sign = util.sign
local lerp = bernstein.lerp1
local lerp2 = bernstein.lerp2
local lerp3 = bernstein.lerp3
local distinct = util.distinct
local significant = util.significant
local negligible = util.negligible
local quadratic = quad.quadratic
local unpack = table.unpack
local getn = table.getn
local pack = table.pack
local floor = math.floor
local ceil = math.ceil
local sqrt  = math.sqrt
local abs  = math.abs
local det2 = util.det2
local atan  = math.atan 
local min = math.min
local max = math.max
local _M = driver.new()
local FLT_MIN = 1.17549435E-38 
local inf   = math.huge

--[[
Funções de Misc(Auxiliares)

]]

local function intersect(x0, y0, x1, y1, x2, y2, x3, y3)
    -- translate x0,y0 to 0,0
    x1, y1 = x1-x0, y1-y0
    x2, y2 = x2-x0, y2-y0
    x3, y3 = x3-x0, y3-y0
    -- find intersection of tangents
    local x10, y10 = x1, y1
    local x32, y32 = x3-x2, y3-y2
    local x30, y30 = x3, y3
    local xp, yp
    -- if tangent at starting point is zero
    if negligible(x10) and negligible(y10) then
        -- if tangent at end point is zero
        if negligible(x32) and negligible(y32) then
            -- cubic degenerates to line
            --error("degenerate to line!")
        else
            -- use other control point of other tangent
            return x2+x0, y2+y0
        end
    -- if tangent at end point zero
    elseif negligible(x32) and negligible(y32) then
        -- use other control point of other tangent
        return x1+x0, y1+y0
    -- actually compute intersection
    else
        local t = det2(x30,x32,y30,y32)/det2(x10,x32,y10,y32)
        return t*x10+x0, t*y10+y0
    end
end

local function crosspmatrix(x0, y0, x1, y1, x2, y2, x3, y3)
    local u0, u1, u2, u3 = power3(x0, x1, x2, x3)
    local v0, v1, v2, v3 = power3(y0, y1, y2, y3)
    local d1 = 0.
    local d2 = det2(u2, u3, v2, v3)
    local d3 = -det2(u1, u3, v1, v3)
    local d4 = det2(u1, u2, v1, v2)
    return d1, d2, d3, d4
end

local function posfunc(n)
	if n > 0 then
		return 1
	elseif n < 0 then 
		return 0
	else
		return 0
	end
end

local function iszero(n)
	if n == 0 then
		return 1
	else
		return 0
	end
end

local function nposfunc(n)
	if n > 0 then
		return 0
	elseif n <= 0 then 
		return 1
	end
end


local function atan2(y,x)
	local temp = 0 
	if x>0 and y>0 then
		temp = atan(y/x)	
	elseif  x<0 and y>0 then
		temp = atan(y/x)
	elseif  x<0 and y<0 then
		temp = atan(y/x)+math.pi
	elseif  x>0 and y<0 then
		temp = atan(y/x)+math.pi
	elseif x == 0 then 
		temp = math.pi*sign(y)/2
	elseif y == 0 then 
		temp = math.pi*sign(x)
	end
	return 180*temp/math.pi
end


function linear_coefficients(x0,x1)
	--print(x0,x1)
	return x1 - x0,x0
end
function quadratic_coefficients(x0,x1,x2)
	return (x0 - 2*x1 + x2), 2*(x1-x0), x0
end
function cubic_coefficients(x0,x1,x2,x3)
	return (x3+3*x1-3*x2-x0), 3*(x0-2*x1+x2), 3*(x1-x0), x0
end

--[[
Funções de Debug

]]



-- these are the two functions that you need to modify/implement
function dump(tbl, indent)
	if not indent then indent = 1 end
	space = "   "
	print(string.rep(space,indent-1) .. "{")
	for k, v in pairs(tbl) do
		formatting = string.rep(space, indent) .. k .. ": "
		if type(v) == "table" then
			print(formatting)
			dump(v, indent+1)
		elseif(type(v) ~= 'function')then
			print(formatting .. v)
		elseif(type(v) == 'image')then
			print(formating ..image.png.string8(v))
		else
			print(formating .. type(v))
		end
	end
	print(string.rep(space,indent-1) .. "}")
end 

--[[
Funções de Converção de Tipo

]]

function makecircle(cx, cy, r,xf)
    -- we start with a unit circle centered at the origin
    -- it is formed by 4 arcs covering each quarter of the unit circle (Using Four Vertex Theorem[Manfredo C.] )
    local s = sqrt(2)/2          -- sin(pi/2)
    local c = sqrt(2)/2 -- cos(pi/2)
    local w = sqrt(2)/2
    local newpath = _M.path{
        _M.M,  1,  0,
        _M.R,  c,  s,  w, 0, 1,
        _M.R,  -c, s,  w,  -1, 0,
        _M.R,  -c, -s,  w,  0,  -1,
	_M.R,  c,  -s,  w,  1,0,		-- colocar em baixo os parametros do circulo
        _M.Z}:scale(r, r):translate(cx, cy)
	newpath.xf = xf*newpath.xf	
    return newpath
end


function maketriangle(x1, y1, x2, y2, x3, y3)
	local newpath = _M.path{_M.M,x1,y1,_M.L,x2,y2,_M.L,x3,y3,_M.L,x1,y1,_M.Z}
	newpath.xmin = min(x1,x2,x3)
	newpath.ymin = min(y1,y2,y3)
	newpath.xmax = max(x1,x2,x3)
	newpath.ymax = max(y1,y2,y3)
	--dump(newpath)
	return newpath
end

function makepolygon(data)
	local newpath = _M.path()
	newpath:open()
	newpath:begin_closed_contour(_,data[1],data[2])
	local px,py = data[1],data[2]
	newpath.xmin = px
	newpath.ymin = py
	newpath.xmax = px
	newpath.ymax = py
	for i = 3, #data, 2 do 
		newpath:linear_segment(px,py,data[i],data[i+1])
		px,py = data[i],data[i+1]
		newpath.xmin = min(newpath.xmin,px)
		newpath.ymin = min(newpath.ymin,py)
		newpath.xmax = max(newpath.xmax,px)	
		newpath.ymax = max(newpath.ymax,py)
	end
	--dump(newpath)
	newpath:linear_segment(px,py,data[1],data[2])
	newpath:end_closed_contour(data[1],data[2])
    return newpath
end



-- cut canonic rational quadratic segment and recanonize
local function cutr2s(a, b, x0, y0, x1, y1, w1, x2, y2)
    local u0 = lerp2(a,a,x0, x1, x2)
    local v0 = lerp2(a,a,y0, y1, y2)
    local r0 = lerp2(a,a,1, w1, 1)
    local u1 = lerp2(a,b,x0, x1, x2)
    local v1 = lerp2(a,b,y0, y1, y2)
    local r1 = lerp2(a,b,1, w1, 1)
    local u2 = lerp2(b,b,x0, x1, x2)
    local v2 = lerp2(b,b,y0, y1, y2)
    local r2 = lerp2(b,b,1, w1, 1)
    local ir0, ir2 = 1/r0, 1/r2
    assert(ir0*ir2 >= 0, "canonization requires split!")
    local ir1 = sqrt(ir0*ir2)
    return u0*ir0, v0*ir0, u1*ir1, v1*ir1, r1*ir1, u2*ir2, v2*ir2
end

-- here are functions to find a root in a rational quadratic
-- you can write your own functions to find roots for lines,
-- integral quadratics, and cubics

-- assumes monotonic and x0 <= 0 <= x2
local function recursivebisectrationalquadratic(x0, x1, w1, x2, ta, tb, n)
    local tm = 0.5*(ta+tb)
    local xm = lerp2( tm, tm,x0, x1, x2)
    local wm = lerp2( tm, tm,1, w1, 1)
    if abs(xm) < TOL*abs(wm) or n >= MAX_ITER then
        return tm
    else
        n = n + 1
        if (wm < 0) ~= (xm < 0) then ta = tm
        else tb = tm end
        -- tail call
        return recursivebisectrationalquadratic(x0, x1, w1, x2, ta, tb, n)
    end
end

-- assumes monotonic and root in [0, 1]
local function bisectrationalquadratic(x0, x1, w1, x2)
    -- ensure root is bracketed by [0,1]
    assert(x2*x0 <= 0, "no root in interval!")
    -- reorder segment so it is increasing
    if x0 > x2 then
        return 1-recursivebisectrationalquadratic(x2, x1, w1, x0, 0, 1, 0)
    else
        return recursivebisectrationalquadratic(x0, x1, w1, x2, 0, 1, 0)
    end
end


local function recursivebisectline(x0, x1, ta, tb, n)
    local tm = 0.5*(ta+tb)
    local xm = lerp(tm,x0, x1)
    if abs(xm) < TOL or n >= MAX_ITER then
        return tm
    else
        n = n + 1
        if xm < 0 then ta = tm
        else tb = tm end
        -- tail call
        return recursivebisectline(x0, x1, ta, tb, n)
    end
end
-- assumes monotonic and root in [0, 1]
local function bisectline( x0 , x1 )
    assert(x1*x0 <= 0, "no root in interval!")
    if x0 > x1 then
        return 1-recursivebisectline(x1, x0, 0, 1, 0)
    else
        return recursivebisectline(x0, x1, 0, 1, 0)
    end
end

local function recursivebisectquadratic(x0, x1, x2, ta, tb, n)
    local tm = 0.5*(ta+tb)
    local xm = lerp2( tm, tm,x0, x1, x2)
    if abs(xm) < TOL or n >= MAX_ITER then
        return tm
    else
        n = n + 1
        if xm < 0 then ta = tm
        else tb = tm end
        -- tail call
        return recursivebisectquadratic(x0, x1, x2, ta, tb, n)
    end
end

-- assumes monotonic and root in [0, 1]
local function bisectquadratic( x0 , x1 , x2 )
    assert(x2*x0 <= 0, "no root in interval!")
    if x0 > x2 then
        return 1-recursivebisectquadratic(x2 , x1, x0, 0, 1, 0)
    else
        return recursivebisectquadratic(x0, x1, x2, 0, 1, 0)
    end
end

local function recursivebisectcubic(x0, x1, x2, x3, ta, tb, n)
    local tm = 0.5*(ta+tb)
    local xm = lerp3( tm, tm, tm,x0, x1, x2, x3 )
    if abs(xm) < TOL or n >= MAX_ITER then
        return tm
    else
        n = n + 1
        if xm < 0 then ta = tm
        else tb = tm end
        -- tail call
        return recursivebisectcubic(x0, x1, x2, x3, ta, tb, n)
    end
end
-- assumes monotonic and root in [0, 1]
local function bisectcubic( x0 , x1 , x2 , x3 )
    assert(x3*x0 <= 0, "no root in interval!")
   if x0 > x3 then
        return 1-recursivebisectcubic(x3, x2, x1, x0, 0, 1, 0)
    else
        return recursivebisectcubic(x0, x1, x2, x3, 0, 1, 0)
    end
end

-- transforms path by xf and ensures it is closed by a final segment
local function newxformer(xf, forward)
    local fx, fy -- first contour cursor
    local px, py -- previous cursor
    local xformer = {}
    function xformer:begin_closed_contour(len, x0, y0)
        fx, fy = xf:apply(x0, y0,1)
        forward:begin_closed_contour(_, fx, fy)
        px, py = fx, fy
    end
    xformer.begin_open_contour = xformer.begin_closed_contour
    function xformer:end_closed_contour(len)
        if px ~= fx or py ~= fy then
            forward:linear_segment(px, py, fx, fy)
        end
        forward:end_closed_contour(_)
    end
    xformer.end_open_contour = xformer.end_closed_contour
    function xformer:linear_segment(x0, y0, x1, y1)
       x1, y1 = xf:apply(x1, y1,1)
       forward:linear_segment(px, py, x1, y1)
       px, py = x1, y1
    end
    function xformer:quadratic_segment(x0, y0, x1, y1, x2, y2)
        x1, y1 = xf:apply(x1, y1,1)
        x2, y2 = xf:apply(x2, y2,1)
        forward:quadratic_segment(px, py, x1, y1, x2, y2)
        px, py = x2, y2
    end
    function xformer:rational_quadratic_segment(x0, y0, x1, y1, w1, x2, y2)
        x1, y1, w1 = xf:apply(x1, y1, w1)
        x2, y2 = xf:apply(x2, y2,1)
        --assert(w1 > FLT_MIN, "unbounded rational quadratic segment")
        forward:rational_quadratic_segment(px, py, x1, y1, w1, x2, y2)
        px, py = x2, y2
    end
    function xformer:cubic_segment(x0, y0, x1, y1, x2, y2, x3, y3)
        x1, y1 = xf:apply(x1, y1,1)
        x2, y2 = xf:apply(x2, y2,1)
        x3, y3 = xf:apply(x3, y3,1)
        forward:cubic_segment(px, py, x1, y1, x2, y2, x3, y3)
        px, py = x3, y3
    end
    function xformer:degenerate_segment(x0, y0, dx0, dy0, dx1, dy1, x1, y1)
      	x1, y1 = xf:apply(x1, y1)
       dx0, dy0 = xf:apply(dx0, dy0, 0)
       dx1, dy1 = xf:apply(dx1, dy1, 0)
       forward:degenerate_segment(px, py, dx0, dy0, dx1, dy1, x1, y1)
       px, py = x1, y1	
    end
    return xformer
end	
-- remove segments that degenerate to points
-- should be improved to remove "almost" points
-- assumes segments are monotonic
local function newcleaner(forward)
    --dump(forward)
    local cleaner = {}
    function cleaner:begin_closed_contour(len, x0, y0)
        forward:begin_closed_contour(_, x0, y0)
    end
    cleaner.begin_open_contour = cleaner.begin_closed_contour
    function cleaner:linear_segment(x0, y0, x1, y1)
        if x0 ~= x1 or y0 ~= y1 then
            forward:linear_segment(px, py, x1, y1)
        end
    end
    function cleaner:quadratic_segment(x0, y0, x1, y1, x2, y2)
        if x0 ~= x2 or y0 ~= y2 then
            forward:quadratic_segment(x0, y0, x1, y1, x2, y2)
        end
    end
    function cleaner:rational_quadratic_segment(x0, y0, x1, y1, w1, x2, y2)
        if x0 ~= x2 or y0 ~= y2 then
            if abs(w1-1) > FLT_MIN then
                forward:rational_quadratic_segment(x0, y0, x1, y1, w1, x2, y2)
            else
                forward:quadratic_segment(x0, y0, x1, y1, x2, y2)
            end
        end
    end
    function cleaner:cubic_segment(x0, y0, x1, y1, x2, y2, x3, y3)
        if x0 ~= x3 or y0 ~= y3 then
            forward:cubic_segment(x0, y0, x1, y1, x2, y2, x3, y3)
        end
    end
    function cleaner:degenerate_segment(x0, y0, dx0, dy0, dx1, dy1, x1, y1)
	forward:degenerate_segment(x0, y0, dx0, dy0, dx1, dy1, x1, y1)	
    end
    function cleaner:end_closed_contour(len)
        forward:end_closed_contour(_)
    end
    cleaner.end_open_contour = cleaner.end_closed_contour
    return cleaner
end	

local function newmonotonizer(forward)
    local monotonizer = {}
    function monotonizer:begin_closed_contour(len, x0, y0)
        forward:begin_closed_contour(_, x0, y0)
    end
    monotonizer.begin_open_contour = monotonizer.begin_closed_contour
    function monotonizer:linear_segment(x0, y0, x1, y1)
        forward:linear_segment(x0, y0, x1, y1) -- o segmento linear não precisa ser monotonizado
    end
    function monotonizer:quadratic_segment(x0, y0, x1, y1, x2, y2)
        --descobre as raízes de x'(t) e y'(t) ordena os t's e usa lerp2 pra descobrir os pontos de controle      
	--print(x0, y0, x1, y1, x2, y2)
	local l1,l2,l3,l4 = x1-x0,y1-y0 ,x2-x0,y2-y0
	local e = l3*l2 - l1*l4
	if e ~= 0 then    
        local t = { 0 } -- valores de t para os pontos que representam os segmentos monotônicos 
    
        if ( x0 + x2 ~= 2*x1 ) then
	    --print("Passei Aqui")
            -- caso a raiz não caia no intervalo [0,1], o resultado não nos interessa
            if ( (x0 - x1)/(x0 - 2*x1 + x2) < 1 and (x0 - x1)/(x0 - 2*x1 + x2) > 0 ) then 
                t[#t + 1] =  (x0 - x1)/(x0 - 2*x1 + x2)--raiz de x'(t) = 0
            end
            
        end
	--print("Passei Aqui")
        if ( y0 + y2 ~= 2*y1 ) then
		--print("Passei Aqui")
            -- caso a raiz não caia no intervalo [0,1], o resultado não nos interessa
            if ( (y0 - y1)/(y0 - 2*y1 + y2) < 1 and (y0 - y1)/(y0 - 2*y1 + y2) > 0 ) then 
                t[#t + 1] =  (y0 - y1)/(y0 - 2*y1 + y2)--raiz de y'(t) = 0
            end
        end
	
        t[#t + 1] = 1
        --coloca os t's em ordem crescente ( Quick Sort)
        table.sort(t, quicksort)
        for i = 1, (#t - 1)  do
            local px0 = lerp2(t[i],t[i],x0,x1,x2)
            local py0 = lerp2(t[i],t[i],y0,y1,y2)

            local px1 = lerp2(t[i],t[i+1],x0,x1,x2)
            local py1 = lerp2(t[i],t[i+1],y0,y1,y2)

            local px2 = lerp2(t[i+1],t[i+1],x0,x1,x2)
            local py2 = lerp2(t[i+1],t[i+1],y0,y1,y2)
	    
            --já pode dar o foward do quadratic segment pra esses 2 junto com o anterior
            forward:quadratic_segment(px0,py0,px1,py1,px2,py2)
        end
	else
	    forward:linear_segment(x0, y0, x2, y2)
	end
   end
    function monotonizer:rational_quadratic_segment(x0, y0, x1, y1, w1, x2, y2)
	local l1,l2,l3,l4 = x1/w1-x0,y1/w1-y0 ,x2-x0,y2-y0
	local e = l3*l2 - l1*l4
	--print(x0, y0, x1, y1, w1, x2, y2,s,a1,e)
	if e ~= 0 then 
		local r = {0,1}
		local a,b,c = quadratic_coefficients(y0,y1,y2)
		local d,e,f = quadratic_coefficients(1,w1,1) 
		local  dca, dcb, dcc = (a*e - b*d), 2*(a*f-c*d), b*f - c*e
		local root = {quadratic(dca,dcb,dcc)}
		if root[1] == 2 then
			for i = 2,4,2 do
				local t = root[i]/root[i+1]
				if t > 0 and t < 1 then
					table.insert(r,t)
				end
			end
		end
		local a,b,c = quadratic_coefficients(x0,x1,x2)
		local d,e,f = quadratic_coefficients(1,w1,1) 
		local  dca, dcb, dcc = (a*e - b*d), 2*(a*f-c*d), b*f - c*e
		local root = {quadratic(dca,dcb,dcc)}
		if root[1] == 2 then
			for i = 2,4,2 do
				local t = root[i]/root[i+1]
				if t > 0 and t < 1 then
					table.insert(r,t)
				end
			end
		end
		table.sort(r,quicksort)
		for i = 1, #r - 1 do
			 forward:rational_quadratic_segment(cutr2s(r[i],r[i+1],x0,y0,x1,y1,w1,x2,y2))
		 end
	else
	    forward:linear_segment(x0, y0, x2, y2)
	end
    end
    function monotonizer:cubic_segment(x0, y0, x1, y1, x2, y2, x3, y3)
	local l1,l2,l3,l4,l5,l6 = x1-x0,y1-y0 ,x2-x0,y2-y0,x3-x0,y3-y0
	local a1 = l4
	local e1 = l3*l2 - l1*l4
	local e2 = l5*l4 - l3*l6
	local e3 = l5*l2 - l1*l6
	local s = sign(-2*l6*e1*e2)
	local d1,d2,d3,d4 = crosspmatrix(x0, y0, x1, y1, x2, y2, x3, y3)
	if d2 ~= 0 or d3 ~=0  then
		-- raciocínio análogo ao quadratic_segment
		local t = { 0 } -- valores de t para os pontos que representam os segmentos monotônicos 
		local Qx = {} -- vetor da coordenada x dos novos pontos de controle
		local Qy = {} -- vetor da coordenada y dos novos pontos de controle
		local solution
		local t1, s1 , t2 , s2
		local r1,r2 = 0,0
		--Teste dos Ponto de Infleção
		--Teste dos Ponto de Duplo
		local d = 3.*d3*d3-4.*d2*d4		
		if d < 0 then
			solution , t1 , s1 , t2 , s2 = quadratic(d2*d2, -d2*d3, d3*d3-d2*d4, -.25*d2*d2*d)
			if ( solution == 2 ) then
			    if ( t1/s1 > 0 and t1/s1 < 1 )  then
				r1 = t1/s1
				r2 = (d3/d2-r1)
			    elseif ( t2/s2 > 0 and t2/s2 < 1 and t2/s2 ~= t1/s1)  then
				 r1 = t2/s2
				 r2 = (d3/d2-r1)
			    end 
			 end
		else
			solution , t1 , s1 , t2 , s2 = quadratic(-3.*d2, 3.*d3, -d4, .25*3.*d)
			if ( solution == 2 ) then
			    if ( t1/s1 > 0 and t1/s1 < 1 )  then
				t[#t + 1] = t1/s1
			    end 
			    if ( t2/s2 > 0 and t2/s2 < 1 and t2/s2 ~= t1/s1)  then
				t[#t + 1] = t2/s2
			    end 
			end
		end
		--print(r1,r2)
		if r1 > 0 and r1 < 1 and r2 > 0 and r2 < 1 then
			t[#t + 1] = r1
			t[#t + 1] = r2
		end
		-- Teste das raizes
		 solution , t1 , s1 , t2 , s2 = quadratic( x3 - 3*x2 + 3*x1 - x0 , 2*x2 - 4*x1 + 2*x0 , x1 - x0 )
		 if ( solution == 2 ) then
		    if ( t1/s1 > 0 and t1/s1 < 1 )  then
		        t[#t + 1] = t1/s1
		    end 
		    if ( t2/s2 > 0 and t2/s2 < 1 and t2/s2 ~= t1/s1)  then
		        t[#t + 1] = t2/s2
		    end 
		 end

		 solution , t1 , s1 , t2 , s2 = quadratic( y3 - 3*y2 + 3*y1 - y0 , 2*y2 - 4*y1 + 2*y0 , y1 - y0 )
		 if ( solution == 2) then
		    if ( t1/s1 > 0 and t1/s1 < 1 )  then
		        t[#t + 1] = t1/s1
		    end 
		    if (  t2/s2 > 0 and t2/s2 < 1 and t2/s2 ~= t1/s1 )  then
		        t[#t + 1] = t2/s2
		    end 
		 end
		 t[#t + 1] = 1
		--dump(t)
		--coloca os t's em ordem crescente ( Quick Sort)
		table.sort(t, quicksort) -- Problema No sort
		--dump(t)
		Qx[1] = x0
		Qy[1] = y0
		for i = 1, #t - 1  do
		    Qx[#Qx + 1] = lerp3(t[i],t[i],t[i+1],x0,x1,x2,x3)
		    Qx[#Qx + 1] = lerp3(t[i],t[i+1],t[i+1],x0,x1,x2,x3)
		    Qx[#Qx + 1] = lerp3(t[i+1],t[i+1],t[i+1],x0,x1,x2,x3)

		    Qy[#Qy + 1] = lerp3(t[i],t[i],t[i+1],y0,y1,y2,y3)
		    Qy[#Qy + 1] = lerp3(t[i],t[i+1],t[i+1],y0,y1,y2,y3)
		    Qy[#Qy + 1] = lerp3(t[i+1],t[i+1],t[i+1],y0,y1,y2,y3)

		    --já pode dar o foward do cubic segment pra esses 3 junto com o anterior
		    forward:cubic_segment(Qx[#Qx - 3], Qy[#Qy - 3], Qx[#Qx - 2], Qy[#Qy - 2], Qx[#Qx - 1], Qy[#Qy - 1], Qx[#Qx], Qy[#Qy])
        end
	elseif d4 ~= 0 then 
		local xi,yi = intersect(x0, y0, x1, y1, x2, y2, x3, y3)
		self.quadratic_segment(self,x0, y0,xi, yi,x3, y3)
	else
		self.linear_segment(self,x0, y0,x3, y3)
	end
    end	
    function monotonizer:end_closed_contour(len)
        forward:end_closed_contour(_)
    end
    function monotonizer:degenerate_segment(x0, y0, dx0, dy0, dx1, dy1, x1, y1)
        --forward:degenerate_segment(x0, y0, dx0, dy0, dx1, dy1, x1, y1)
    end
    monotonizer.end_open_contour = monotonizer.end_closed_contour
    return monotonizer
end

function signalcubic(l1,l2,l3,l4,l5,l6)
	local a = (l2-l4-l6)
	local b = (4*l2^2-2*l2*l4+l4^2)*l5^2
	local c = (9*l4^2-6*l4*l6-4*l6^2)*l1^2
	local d = (9*l2^2-12*l2*l6-l6^2)*l3^2
	local e = -2*l3*(l5*(3*l2^2-l4*l6+l2*(-6*l4+l6))+l1*(l2*(9*l4-3*l6)-l6*(6*l4+l6)))
	local f = 2*l1*l5*(-l4*(6*l4+l6)+l2*(3*l4+4*l6))
	local result =	sign(a*(-b+c+d+f+e))
	return result
end

function transformpath(oldpath,k,elements)
    local xmin,xmax,ymin,ymax = inf,0,inf,0
    local newpath = _M.path()
    local px,py
    newpath:open()
    oldpath:iterate(
        newxformer(oldpath.xf,
          newmonotonizer(
                newcleaner(
                    newpath)))  )	
    newpath:close()
    	newpath.parameter = {}
	-- Gerar Bound Box e os parametros 
	for i,v in ipairs(newpath.instructions) do
	  local o = newpath.offsets[i]
	  local data = newpath.data			  	
	  	  if v == "begin_closed_contour" then
			  px,py = data[o+1],data[o+2]
		  	  newpath.parameter[i] = {0,0}
			  xmin = min(xmin,data[o+1])
		  	  ymin = min(ymin,data[o+2])
			  xmax = max(xmax,data[o+1])
			  ymax = max(ymax,data[o+2])
		  elseif v == "linear_segment" then
			  local a = data[o+3]-data[o+1]
			  local b = data[o]-data[o+2]
			  local c = -a*data[o]-b*data[o+1]
			  newpath.parameter[i] = {a,b,c}
			  xmin = min(xmin,data[o+2],data[o])
			  ymin = min(ymin,data[o+3],data[o+1])
			  xmax = max(xmax,data[o+2],data[o])
			  ymax = max(ymax,data[o+3],data[o+1])
		  elseif v == "begin_open_contour" then
			  px,py = data[o+1],data[o+2]
			  newpath.parameter[i] = {0,0}
			  xmin = min(xmin,data[o+1])
			  ymin = min(ymin,data[o+2])
			  xmax = max(xmax,data[o+1])
			  ymax = max(ymax,data[o+2])	
		  elseif v == "quadratic_segment" then
			x1,y1,x2,y2,x3,y3 = data[o],data[o+1],data[o+2],data[o+3],data[o+4],data[o+5]
			local l1,l2,l3,l4 = x2-x1,y2-y1 ,x3-x1,y3-y1
			local a1 = l4
			local b1 = -l3
			local a2 = l4-l2
			local b2 = l1-l3
			local a = 2*l1 - l3
			local b = l4-2*l2
			local c = l2
			local d = -l1
			local e = l3*l2 - l1*l4
			local s = sign(-2*a1*e)
			newpath.parameter[i] = {s,a1,b1,a2,b2,a,b,c,d,e}
			  xmin = min(xmin,data[o+2],data[o+4],data[o])
			  ymin = min(ymin,data[o+3],data[o+5],data[o+1])
			  xmax = max(xmax,data[o+2],data[o+4],data[o])
			  ymax = max(ymax,data[o+3],data[o+5],data[o+1])
		  elseif v == "cubic_segment" then
			  x1,y1,x2,y2,x3,y3,x4,y4 = data[o],data[o+1],data[o+2],data[o+3],data[o+4],data[o+5],data[o+6],data[o+7]
			  local l1,l2,l3,l4,l5,l6 = x2-x1,y2-y1 ,x3-x1,y3-y1,x4-x1,y4-y1
			  -- Calculo do Sinal
			  local la0,lb0 = l2,-l1
			  local la1,lb1,lc1 = l4-l2,l1-l3,l3*l2-l4*l1 
			  local la2,lb2,lc2 = l6-l4,l3-l5,l5*l4-l6*l3
			  local la3,lb3 = l6,-l5
			  local a = (-9*l3*l2+3*l5*l2+9*l1*l4-3*l1*l6)
			  local b1 = (6*l1-3*l3)
			  local b2 = (-6*l2+3*l4)
			  local c1 = (-3*l1+3*l3-l5)
			  local c2 = (3*l2-3*l4+l6)
			  local d = (9*l3*l2-6*l5*l2-9*l1*l4+3*l5*l4+6*l1*l6-3*l3*l6)
		 	  local e1 = -3*l1
			  local e2 = 3*l2
			  local f1 = (-3*l1+3*l3-l5) 
			  local f2 = 9*l3*l2-9*l1*l4 
			  local f3 = (3*l2-3*l4+l6)
			  local s = signalcubic(l1,l2,l3,l4,l5,l6)	 
			  newpath.parameter[i] = {s,la0,lb0,la1,lb1,lc1,la2,lb2,lc2,la3,lb3,a,b1,b2,c1,c2,d,e1,e2,f1,f2,f3}
			  xmin = min(xmin,data[o+2],data[o+4],data[o+6],data[o])
			  ymin = min(ymin,data[o+3],data[o+5],data[o+7],data[o+1])
			  xmax = max(xmax,data[o+2],data[o+4],data[o+6],data[o])
			  ymax = max(ymax,data[o+3],data[o+5],data[o+7],data[o+1])
		  elseif v == "rational_quadratic_segment" then
			  x1,y1,x2,y2,w2,x3,y3 = data[o],data[o+1],data[o+2],data[o+3],data[o+4],data[o+5],data[o+6]
			  local l1,l2,l3,l4 = (x2-x1*w2),(y2-y1*w2),(x3-x1),(y3-y1)
		          local s =  (1+w2)*l4*(l1*l4-l3*l2)
			  local a0 = l2/w2
			  local b0 = -l1/w2
			  local a1 = l4
			  local b1 = -l3
		 	  local a2 = l4-l2/w2
			  local b2 = l1/w2-l3
			  local c2 = l3*l2/w2 - l1*l4/w2
			  local a = (4*l1^2-4*w2*l1*l3+l3^2)
			  local b = 4*l1*l3*l2-4*l4*l1^2	
			  local c = (-8*l1*l2+4*w2*(l3*l2+l1*l4)-2*l3*l4)
			  local d = (4*l2^2-4*w2*l2*l4+l4^2)
			  local e = -4*l3*l2^2+4*l1*l2*l4
			  newpath.parameter[i] = {s,a0,b0,a1,b1,a2,b2,c2,a,b,c,d,e,f}
			  xmin = min(xmin,data[o+2],data[o+5],data[o])
			  ymin = min(ymin,data[o+3],data[o+6],data[o+1])
			  xmax = max(xmax,2*data[o+2]/(1+w2),data[o+5],data[o])
			  ymax = max(ymax,2*data[o+3]/(1+w2),data[o+6],data[o+1])
		  elseif v == "end_open_contour" then
			  local a = data[o+3]-data[o+1]
			  local b = data[o]-data[o+2]
			  local c = -a*data[o]-b*data[o+1]
			  newpath.parameter[i] = {a,b,c}
			  xmin = min(xmin,data[o])
			  ymin = min(ymin,data[o+1])
			  xmax = max(xmax,data[o])
			  ymax = max(ymax,data[o+1])
		  elseif v == "end_closed_contour" then
			  local a = py-data[o+1]
			  local b = data[o]-px
			  local c = -a*data[o]-b*data[o+1]
			  newpath.parameter[i] = {a,b,c}
			  xmin = min(xmin,data[o])
			  ymin = min(ymin,data[o+1])
			  xmax = max(xmax,data[o])
			  ymax = max(ymax,data[o+1])
		  elseif v == "degenerate_segment" then
			  xmin = min(xmin,data[o+2],data[o+4],data[o])
			  ymin = min(ymin,data[o+3],data[o+5],data[o+1])
			  xmax = max(xmax,data[o+2],data[o+4],data[o])
			  ymax = max(ymax,data[o+3],data[o+5],data[o+1])	
		end
	end
	if oldpath.xmin == nil and oldpath.xmax == nil and oldpath.ymin == nil and oldpath.ymax == nil then 
		newpath.xmin = xmin
		newpath.xmax = xmax 
		newpath.ymin = ymin
		newpath.ymax = ymax
    	else 
		newpath.xmin = oldpath.xmin
		newpath.xmax = oldpath.xmax 
		newpath.ymin = oldpath.ymin
		newpath.ymax = oldpath.ymax
	end	
    return newpath
end
--[[
Funções de Elipses

]]
local function applytheconic(conic,x,y)
	local v = vector.vector(x,y)
	local aux =  conic * v
	local result,rest = v:dot(aux)
	return result + rest
end

local function conic2origin(conic)
	local a,h,g,h1,b,f,g1,f1,c = unpack(conic)
	local p = (-b*g + f*h)/(a*b - h*h)
	local q = (a*f - g*h)/(-a*b + h*h)
	return p,q
end

local function boundBox(shape, conic)
	-- translate the center of ellipse to the origin
	t = xform.translate(conic2origin(conic))
	local nconic = t:transpose()*conic*t
	-- new coefficients
	a,h,g,h1,b,f,g1,f1,c = unpack(nconic)
	-- limits of the box
	y1 = sqrt((a*c)/ (h*h - a*b))
	y2 = -1*y1

	x1 = sqrt((b*c)/ (h*h - a*b))
	x2 = -1*x1

	-- needs the inverse transfomation
	x1,y1 = t:apply(x1,y1)
	x2,y2 = t:apply(x2,y2)
	return x2,y2,x1,y1 
end

--[[
Funções de Preparo de Cena

]]

local function preparescene(scene)
	local xmin,ymin,xmax,ymax = inf,inf,0,0
   	local xf = scene.xf
	for k = #scene.elements,1,-1 do
		e = scene.elements[k]
		e.shape.xf = xf * e.shape.xf
		e.paint.xf = xf * e.paint.xf
		if e.shape['type'] == 'circle' then
			local s = e.shape
			local a,b,f,g = 1,1,-1*s.cy,-1*s.cx
			local c = s.cx*s.cx + s.cy*s.cy - s.r*s.r
			s.conic = xform.xform(a,0,g, 0,b,f, g,f,c)
			s.conic = s.xf:inverse():transpose() * s.conic * s.xf:inverse() 
			xmin,ymin,xmax,ymax = boundBox(s,s.conic)
			e.shape = makecircle(s.cx, s.cy, s.r,e.shape.xf)
			e.shape.xmin = xmin
			e.shape.ymin = ymin
			e.shape.xmax = xmax
			e.shape.ymax = ymax
		elseif e.shape['type'] == 'triangle' then
			local s = e.shape
			e.shape.x1,e.shape.y1 =e.shape.xf:apply(e.shape.x1,e.shape.y1)
			e.shape.x2,e.shape.y2 =e.shape.xf:apply(e.shape.x2,e.shape.y2)
			e.shape.x3,e.shape.y3 =e.shape.xf:apply(e.shape.x3,e.shape.y3)
			e.shape = maketriangle(s.x1, s.y1, s.x2, s.y2, s.x3, s.y3)
		elseif e.shape['type'] == 'polygon' then
			for i = 1,#e.shape.data,2 do
				e.shape.data[i],e.shape.data[i+1]=e.shape.xf:apply(e.shape.data[i],e.shape.data[i+1])
				xmin = min(xmin,e.shape.data[i])
				ymin = min(ymin,e.shape.data[i+1])
				xmax = max(xmax,e.shape.data[i])
				ymax = max(ymax,e.shape.data[i+1])
			end
			e.shape = makepolygon(e.shape.data)
			e.shape.xmin = xmin
			e.shape.ymin = ymin
			e.shape.xmax = xmax
			e.shape.ymax = ymax
		end
		if e.shape['type'] == 'path' then
			e.shape = transformpath(e.shape,k,scene.elements)
		end
		--dump(e)
		if e.paint['type'] == 'texture' then
			e.paint.xf = e.paint.xf:inverse()
		end
		if e.paint['type'] == 'lineargradient' then
			local p = e.paint

			local tp1 = xform.translate(unpack(p.data.p1)):inverse()
			--p.data.tp1 = tp1
			p.data.tp2 = {tp1:apply(unpack(p.data.p2))}

			local degree = atan2(p.data.tp2[2],p.data.tp2[1])
			local rot = xform.rotate(-degree)

			-- rotate p2 to be in the x-axis
			p.data.tp2 = {rot:apply(tp1:apply(unpack(p.data.p2)))}
			p.data.px = p.data.tp2[1]
			p.xf = rot*tp1*p.xf:inverse()
		end

		if e.paint['type'] == 'radialgradient' then
			local p = e.paint

			local center = p.data.center
			local r = p.data.radius

			-- use implicity representation
			local a,b,f,g = 1,1,-center[2],-center[1]
			local c = center[1]*center[1] + center[2]*center[2] - r*r
			p.circle = xform.xform(a,0,g, 0,b,f, g,f,c)

			-- translate the focus to the origin
			local tfocus = xform.translate(unpack(p.data.focus)):inverse()

			-- translate the focus to the origin, center and the circle
			p.data.tcenter = {tfocus:apply(unpack(p.data.center))}
			p.circle = tfocus:inverse():transpose() * p.circle * tfocus:inverse() 
			--dump(e)
			assert(applytheconic(p.circle,p.data.tcenter[1],p.data.tcenter[2]) < 0, 
			"the center is out of the conic")
			--dump(p.data.tcenter)
			if not(p.data.tcenter[2] == 0 and p.data.tcenter[1] == 0) then
				local degree = atan2(p.data.tcenter[2],p.data.tcenter[1])
				local rot = xform.rotate(-degree)

				p.data.tcenter = {rot:apply(unpack(p.data.tcenter))}
				p.circle = rot:inverse():transpose() * p.circle * rot:inverse() 

				local centerscale = 1/p.data.tcenter[1]
				local tscale = xform.scale(centerscale)

				p.data.tcenter = {tscale:apply(unpack(p.data.tcenter))}
				p.circle = tscale:inverse():transpose() * p.circle * tscale:inverse()

				
				p.circleRadius = sqrt(abs(p.circle[3+6]/p.circle[1] - 1))
				p.xf = rot * tscale * tfocus * p.xf:inverse()
			else 
				p.circleRadius = sqrt(abs(p.circle[3+6]/p.circle[1] - 1))
				p.xf =   tfocus *p.xf:inverse()			
			end
		end
	end
	return scene
end
--[[
Funções de escolha de cor
]]
local function ajusttoramp(ramp,x)
	local xmin, xmax,xmint, xmaxt,index = 0,0,0,0,0,0
	local r,g,b,a =0,0,0,0
	local result = 0
	for i = 1,#ramp-2,2 do
			xmint = min(xmint,ramp[i])
			xmaxt = max(xmaxt,ramp[i+2])
	if (ramp[i] <= x) and (x <= ramp[i+2]) then
		xmin, xmax = ramp[i], ramp[i+2]
		index = i
	end
	end
	if ramp.spread == "pad" then 				
		if xmax ~= xmin then
			result2 = (x - xmin)/(xmax - xmin)
		elseif x <= xmaxt then
			result2 = 0
		elseif x > xmaxt then
			result2 = 1
		end		
		result = min(1,max(0,result2))
		if result == 0 then	
			index = 1			
		elseif result == 1 then
			index = #ramp-3
		end 		
	elseif ramp.spread == "repeat" then 
		result =  abs(x -floor(x))
		for i = 1,#ramp-2,2 do
		if (ramp[i] <= result) and (result <= ramp[i+2]) then
			xmin, xmax = ramp[i], ramp[i+2]
			index = i
		end
		end
		if result == 0 then	
			index = 0			
		elseif result == 1 then
			index = 0
		end 				
	elseif ramp.spread == "reflect" then 
		result =  x -floor(x)
		for i = 1,#ramp-2,2 do
		if (ramp[i] <= result) and (result <= ramp[i+2]) then
			xmin, xmax = ramp[i], ramp[i+2]
			index = i
		end
		end		
		result2 = (x - xmin)/(xmax - xmin)
		result = 2*abs(result2/2 - (floor((result2+1)/2)))
		if result == 0 then	
			index = 0			
		elseif result == 1 then
			index = 0
		end 
	elseif ramp.spread == "transparent" then 
		if xmax ~= xmin then
			result2 = (x - xmin)/(xmax - xmin)
		elseif x <= xmaxt then
			result2 = 0
		elseif x > xmaxt then
			result2 = 1
		end		
		result = min(1,max(0,result2))
		if result == 0 then	
			index = 0			
		elseif result == 1 then
			index = 0
		end 		
	end
	if index ~= 0 then	
		r = (1 - result)*ramp[index+1][1] + result*ramp[index+3][1]
		g = (1 - result)*ramp[index+1][2] + result*ramp[index+3][2]
		b = (1 - result)*ramp[index+1][3] + result*ramp[index+3][3]
		a = (1 - result)*ramp[index+1][4] + result*ramp[index+3][4]
	end 	
        return r,g,b,a
end

local function ajusttorampsimple(x,xmin,xmax,spread)
	local result2 = 0	
	local result = 0
	if spread == "pad" then 				
		if xmax ~= xmin then
			result2 = (x - xmin)/(xmax - xmin)
		elseif x <= xmaxt then
			result2 = 0
		elseif x > xmaxt then
			result2 = 1
		end		
		result = min(1,max(0,result2))
		if result == 0 then	
			index = 1			
		elseif result == 1 then
			index = #ramp-3
		end 		
	elseif spread == "repeat" then 
		result =  abs(x -floor(x))
		if result == 0 then	
			index = 0			
		elseif result == 1 then
			index = 0
		end 				
	elseif spread == "reflect" then 
		result =  x -floor(x)	
		result2 = (x - xmin)/(xmax - xmin)
		result = 2*abs(result2/2 - (floor((result2+1)/2)))
		if result == 0 then	
			index = 0			
		elseif result == 1 then
			index = 0
		end 
	elseif spread == "transparent" then 
		if xmax ~= xmin then
			result2 = (x - xmin)/(xmax - xmin)
		elseif x <= xmaxt then
			result2 = 0
		elseif x > xmaxt then
			result2 = 1
		end		
		result = min(1,max(0,result2))
		if result == 0 then	
			index = 0			
		elseif result == 1 then
			index = 0
		end 		
	end
	return result
end

function lineargradient(paint,x,y)
	local ramp  = paint.data.ramp
	local opacity = paint.opacity
	local x,y,w = paint.xf:apply(x,y)
	local result  = x/paint.data.px
	local r,g,b,a = ajusttoramp(ramp,result)
	return {r,g,b,	opacity*a}
end

function texturegradient(paint,x,y)
	local r,g,b,a =	0,0,0,0
	local w,h = paint.data.image.width,paint.data.image.height
	local opacity = paint.opacity
	local tx,ty =  paint.xf:apply(x,y)
	local tx,ty = ajusttorampsimple(tx,0,1,paint.data.spread),ajusttorampsimple(ty,0,1,paint.data.spread)
	if tx <=1 and ty <= 1 and tx > 0 and ty > 0 then
		r,g,b,a =  paint.data.image:get(ceil(w*tx),ceil(h*ty))
	end	
	return {r,g,b,opacity*a}
end

function opacity(paint)
	local r,g,b,a = unpack(paint.data)
	local opacity = paint.opacity
	return {r,g,b,	opacity*a}
end

function radialgradient(paint,xi,yi)
	local i = 0
	local x,y,w = paint.xf:apply(xi,yi)
	local r,g,b,a = 0,0,0,0
	local opacity = paint.opacity
	local ramp    = paint.data.ramp
	local t = 0
	if not(paint.data.tcenter[2] == 0 and paint.data.tcenter[1] == 0) then
		local a = x*x + y*y
		local b = -2*x
		local c = 1 - paint.circleRadius*paint.circleRadius
		local root = {quadratic(a,b,c)}
		t = root[3]/root[2]
	elseif  (paint.data.tcenter[2] == 0 and paint.data.tcenter[1] == 0) then
		t = sqrt((x*x + y*y))/(paint.circleRadius)
	end 
	r,g,b,a = ajusttoramp(ramp,t)	
	return {r,g,b,	opacity*a}
end

local function composecolor(ncolor,color,nsamples)
	local r, g, b, a = unpack(ncolor)
	local nr, ng, nb, alpha = unpack(color)
	nr, ng, nb, alpha =nr, ng, nb, alpha/nsamples
	local af =0
	if(a > 0)then
		af = (alpha+a*(1-alpha))
		r = (alpha*nr + (1-alpha)*r*a)/af
		g = (alpha*ng + (1-alpha)*g*a)/af
		b = (alpha*nb + (1-alpha)*b*a)/af
		a = af
	else
		r = nr
		g = ng
		b = nb
		a = alpha
	end	

	return r,g,b,a
end

-- sample scene at x,y and return r,g,b,a
local function sample(scene, x, y)
	local r ,g,b,a = 0,0,0,0
	local winding = 0
	local bg = 0
	local total = #scene.elements
	local shape = 0
	local paint = 0
	local eofill
	local nsamples = #samples/2
	local xmin,xmax,ymin,ymax 
	for j=1,nsamples,1 do
		x = x+samples[2*j-1]
		y = y+samples[2*j]	
		for i= total,1,-1 do
			xmin = scene.elements[i].shape.xmin 
			ymin =  scene.elements[i].shape.ymin 
			xmax =  scene.elements[i].shape.xmax 
			ymax =  scene.elements[i].shape.ymax 
			winding = 0
			if (x-xmin) > 0 and (xmax-x) >= 0 and (y-ymin) > 0 and (ymax-y) >= 0  then 
				shape = scene.elements[i].shape
				paint = scene.elements[i].paint
				eofill = scene.elements[i].type
				if shape.type == "circle" then
					winding = windingcircle(shape,x,y) 
				elseif shape.type == "triangle" then
					winding = windingtriangle(shape,x,y) 
				elseif shape.type == "polygon"  then
					winding = windingpolygon(shape,x,y) 
				elseif shape.type == "path"  then
					winding = windingpath(shape,x,y)
				end
				winding =  math.abs(winding)
				if eofill == "eofill" then
					winding =  winding%2 
				end
				if winding > 0   and paint.type == "solid" then
					r,g,b,a = composecolor(opacity(paint),{r,g,b,a},nsamples)
				elseif winding  > 0 and paint.type == "lineargradient"  then
					r,g,b,a = composecolor(lineargradient(paint,x,y),{r,g,b,a},nsamples)
				elseif winding  > 0 and paint.type == "radialgradient"  then
					r,g,b,a = composecolor(radialgradient(paint,x,y),{r,g,b,a},nsamples)
				elseif winding  > 0 and paint.type == "texture"  then
					r,g,b,a = composecolor(texturegradient(paint,x,y),{r,g,b,a},nsamples)
				end
			end
			if winding  >0 and a == 1 then
				break
			end
		end
	end
	if a < 1 then
		r,g,b,a = composecolor({1,1,1,1},{r,g,b,a},1)
	end
	return  r,g,b,a
end


-- verifies that there is nothing unsupported in the scene
local function checkscene(scene)
    for i, element in ipairs(scene.elements) do
        assert(element.type == "fill" or
               element.type == "eofill", "unsupported element")
        assert(element.shape.type == "circle" or
               element.shape.type == "triangle" or
              element.shape.type == "polygon" or
			   element.shape.type == "path" , "unsuported primitive")
       assert(not element.shape.style, "stroking not unsuported")
       assert(element.paint.type == "solid" or element.paint.type == "texture" or
               element.paint.type == "lineargradient" or
               element.paint.type == "radialgradient", "unsupported paint")
   end
end


--[[
Funções de Checagem Implicita

]]

function implicitline(a,b,c,y1,y2,x,y)
	local i2 = 0
	local s = sign(a)
	a = s*a
	b = s*b
	c = s*c
	if y <= math.max(y1,y2) and y > math.min(y1,y2) and b*(y)+a*(x)+c < 0 then
		i2=  i2+s
	end
	return i2
end

function implicitQuadratic(param,xt,yt)
	local s,a1,b1,a2,b2,a,b,c,d,e = unpack(param)
	local i,i2,i4= 0,0,0,0,0
	local result = 0 
	i2 = implicitline(c,d,0,0,c,xt,yt)+implicitline(a2,b2,e,c,a1,xt,yt)
	i = implicitline(-a1,-b1,0,0,a1,xt,yt) 
	if i+i2 ~= 0  then
		i4 =  s*((a*yt +xt*b)^2 - 4*(xt*c+d*yt)*e)  
		result = sign(a1)*nposfunc(i4)
	elseif  i ~= 0 then
		result = sign(a1) 
	end
	return result
end

function implicitRQuadratic(param,xt,yt)
	local s,a0,b0,a1,b1,a2,b2,c2,a,b,c,d,e,f= unpack(param)
	local i,i2,i4= 0,0,0,0,0
	local result = 0 
	i2 = implicitline(a0,b0,0,0,a0,xt,yt)+implicitline(a2,b2,c2,a0,a1,xt,yt)
	i = implicitline(-a1,-b1,0,0,a1,xt,yt) 
	if i+i2 ~= 0  then
		i4 = (s)*(yt*(a*yt+b)+xt*(e+yt*c+xt*d)) 
		result = sign(a1)*nposfunc(i4)
	elseif  i ~= 0 then
		result = sign(a1) 
	end
	return result
end

function implicitCubic(param,xt,yt)
	local s,la0,lb0,la1,lb1,lc1,la2,lb2,lc2,la3,lb3,a,b1,b2,c1,c2,d,e1,e2,f1,f2,f3 = unpack(param)
	local result =0
	local i,i2,i4= 0,0,0,0,0
	i2 = implicitline(la0,lb0,0,0,la0,xt,yt)+implicitline(la1,lb1,lc1,la0,la1+la0,xt,yt)+implicitline(la2,lb2,lc2,la1+la0,la3,xt,yt)
	i = implicitline(-la3,-lb3,0,0,la3,xt,yt)
	if i+i2 ~= 0  then
		i4 =  -s*(implicitfuncCubic(a,b1,b2,c1,c2,d,e1,e2,f1,f2,f3,xt,yt))  
		result = sign(la3)*nposfunc(i4)
	elseif  i ~= 0 then
		result = sign(la3) 
	end
	return result
end

function implicitfuncCubic(a,b1,b2,c1,c2,d,e1,e2,f1,f2,f3,x,y)
	local b = (b1*y+x*b2)
	local c = (c1*y+x*c2)	
	local e = (e1*y+e2*x)
	local f =  (f1*y+f2+x*f3) 	
	local result = -a*(e*a-b*c)+d*(-b^2+e*f)+c*(b*a-c*f)
	return result
end 					       

--[[
Funções de Contagem de Voltas

]]
function windingpath(shape,x,y)
	local x0,y0x1,y1,x2,y2,x3,y3
	local r ,g,b,a = 0,0,0,0 
	local i2 = 0
	local previous
	--i = #shape.instructions,1,-1 do	
	for i = #shape.instructions,1,-1 do	
		local v = shape.instructions[i]
	  	local o = shape.offsets[i]
		local param = shape.parameter[i]
		if (shape.xmin < x) and (shape.xmax > x) and (shape.ymin < y) and (shape.ymax >= y) then
			if v == "begin_closed_contour" then
			 	x0 = shape.data[o+1]
			 	y0 = shape.data[o+2]
			elseif v == "begin_open_contour" then
			 	x0 = shape.data[o+1]
			 	y0 = shape.data[o+2]
			elseif v == "linear_segment" then
				local y1,y2 = shape.data[o+1],shape.data[o+3]
				i2 = i2 + implicitline(param[1],param[2],param[3],y1,y2,x,y)
			elseif v == "quadratic_segment" then
				x1,y1 = shape.data[o],shape.data[o+1]
				i2 = i2 + implicitQuadratic(param,x-x1,y-y1)
			elseif v == "cubic_segment" then
				x1,y1,x2,y2,x3,y3,x4,y4 = shape.data[o],shape.data[o+1],shape.data[o+2],shape.data[o+3],shape.data[o+4],shape.data[o+5],shape.data[o+6],shape.data[o+7]
				i2 = i2 + implicitCubic(param,x-x1,y-y1)
			elseif v == "rational_quadratic_segment" then
				x1,y1 = shape.data[o],shape.data[o+1]
				i2 = i2 + implicitRQuadratic(param,(x-x1),(y-y1))
			elseif v == "end_open_contour" then
				local y1 = shape.data[o+1]
				i2 = i2 + implicitline(param[1],param[2],param[3],y0,y2,x,y)
			elseif v == "end_close_contour" then
				local y1 = shape.data[o+1]
				i2 = i2 + implicitline(param[1],param[2],param[3],y0,y2,x,y)  
			elseif v == "degenerate_segment" then
				
			  end
			 end
    end
	return i2
end

--[[
Funções Main

]]
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
        for j = 1, width do
            img:set(j, i, sample(scene, j-0.5, i-0.5))
        end
	--	io.write("\n")
--stderr("\n")
    end
stderr("\n")
stderr("rendering in %.3fs\n", time:elapsed())
time:reset()
    -- store output image
    image.png.store8(file, img)
stderr("saved in %.3fs\n", time:elapsed())
end

return _M
