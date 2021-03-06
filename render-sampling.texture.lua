local driver = require"driver"
local image = require"image"
local chronos = require"chronos"

local unpack = table.unpack
local bezier = require("bezier")
local util = require"util"
local blue = require"blue"
local significant = util.significant
local floor, ceil, min, max, sqrt, inf = math.floor, math.ceil, math.min, math.max, math.sqrt, math.huge

local _M = driver.new()
local prepare = {}

---------------------- shortcut operations--------------------------------------
local function sign(x) return (x<0 and -1) or 1 end

function table.empty (self)
    for _, _ in pairs(self) do  return false end
    return true
end

local function truncate_parameter(t)
    if t < 0 or t == -math.huge then t = 0
    elseif t > 1 or t == math.huge then t = 1
    end
    return t
end

local function prod_biP(p,q) -- This is simple product (ax + by + c)*(dx + ey + f)
  local c,b,a = unpack(p)
  local f,e,d = unpack(q)
  return {c*f, a*f + c*d, b*f + c*e, a*e + b*d, a*d, b*e}
end

local function prod_biP_2(p,q) -- Product (ax + by + c)*(dx² + ey² + fxy + gx + hy + i)
  local c,b,a = unpack(p)
  local i,g,h,f,d,e = unpack(q)
  return {c*i, a*i + c*g, b*i + c*h, a*h + b*g + c*f, a*g + c*d, b*h + c*e, a*f + b*d, b*f + a*e, a*d, b*e}
end

local function colinear(x, y)
  if (not x[4]) then
    return (x[1]*(y[2]- y[3]) + x[2]*(y[3] - y[1]) + x[3]*(y[1] - y[2]) == 0)
  else
    return (x[1]*(y[2]- y[3]) + x[2]*(y[3] - y[1]) + x[3]*(y[1] - y[2]) == 0) and
           (x[4]*(y[2]- y[3]) + x[2]*(y[3] - y[4]) + x[3]*(y[4] - y[2]) == 0)
  end
end

local function cross_with_e1(v1, v2, v3) -- cross product with e1 = (1,0,0)
  local a, b, c, d = v1[1], v1[2], v1[3], v1[4]
  local l, e, f, g = v2[1], v2[2], v2[3], v2[4]
  return c*g - d*f, d*e - b*g, b*f - c*e
end

local function mult(s,u)
  local w = {}
  for i=1,#u do w[i] = s*u[i] end
  return w
end

local function addM(u,v)
  assert(#u == #v, "dot product failed")
  local w = {}
  for i=1,#u do w[i] = u[i]+v[i] end
  return w
end

---------------------------Auxiliarie functions---------------------------------
local function dataTransform(shape,data,pos,m) -- I use this to transform and add to the bounding box m points
  for i=0,(m-1) do
    data[pos+2*i], data[pos+(2*i+1)] = shape.xf : apply(data[pos+2*i],data[pos+(2*i+1)])
    shape.xmin = min(shape.xmin or  inf, data[pos+2*i])
  	shape.ymin = min(shape.ymin or  inf, data[pos+2*i+1])
  	shape.xmax = max(shape.xmax or -inf, data[pos+2*i])
  	shape.ymax = max(shape.ymax or -inf, data[pos+2*i+1])
  end
  return data
end

local function blend(f,b) -- color composition
  return {f[4]*f[1] + (1-f[4])*b[4]*b[1],
          f[4]*f[2] + (1-f[4])*b[4]*b[2],
          f[4]*f[3] + (1-f[4])*b[4]*b[3],
               f[4] + (1-f[4])*b[4]      }
end

local function bezier_to_poly2(x,y,w) -- Find the coefficients of the bezier polynomial for rational quadratics
  local w = w or 1
  return {x[1], 2*(x[2]-x[1]), x[1] - 2*x[2] + x[3]},
         {y[1], 2*(y[2]-y[1]), y[1] - 2*y[2] + y[3]},
         {   1,       2*w - 2,              2 - 2*w}
end

local function bezier_to_poly3(x,y) -- Find the coefficients of the bezier polynomial for integral cubics
  return {x[1], -3*x[1] + 3*x[2], 3*x[1] - 6*x[2] + 3*x[3], -x[1] + 3*x[2] - 3*x[3] + x[4]},
         {y[1], -3*y[1] + 3*y[2], 3*y[1] - 6*y[2] + 3*y[3], -y[1] + 3*y[2] - 3*y[3] + y[4]}
end

local function difP(p) -- Return the coefficients of the derivative of a polynomial
  local dp = {0}
  for i=1,(#p-1) do dp[i] = i*p[i+1] end
  return dp
end

local function multP(p,q) -- Return the coefficients of the multiplication of two polynomials
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
    local a = 2*(w - 1)*(x0 - x2)
    local b = 2*(x0 - 2*w*x0 + 2*x1 - x2)
    local c = 2*(w*x0 - x1)

    n, r1, s1, r2, s2 = _M.quadratic(a, b, c)

    local out1, out2 = 0, 0
    if n > 0 then out1 = r1/s1 end
    if n > 1 then out2 = r2/s2 end

    return truncate_parameter(out1), truncate_parameter(out2)
end

local bezier_cut = {0,
  function(a,b,x,y)
    u0, v0, u1, v1, u2, v2 = bezier.cut2(a, b, x[1], y[1], x[2], y[2], x[3], y[3])
    return {u0,u1,u2}, {v0,v1,v2}
  end,
  function(a,b,x,y)
    u0, v0, u1, v1, u2, v2, u3, v3 = bezier.cut3(a, b, x[1], y[1], x[2], y[2], x[3], y[3], x[4], y[4])
    return {u0,u1,u2,u3}, {v0,v1,v2,v3}
  end,
  r = function(a,b,x,y,w)
    u0, v0, u1, v1, r, u2, v2 = bezier.cut2rc(a, b, x[1], y[1], x[2], y[2], w, x[3], y[3])
    return {u0,u1,u2}, {v0,v1,v2}, r
  end
}

local function polyToFun(p)
  local n = (#p-1)
  return function (t)
    local v = p[n+1]
    for i=n,1,-1 do v = t*v + p[i] end
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

local function poly_num_root(p,a,b) -- This return the set of roots in [a,b] of a polynomial with any degree
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
    if p[2] then return {-p[1]/p[2]} end
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

local partition = {0,
  function(x,y)
      local partition = {0,1}
      partition[3] = truncate_parameter( (x[1]-x[2]) / (x[1] - 2*x[2] + x[3]) )
      partition[4] = truncate_parameter( (y[1]-y[2]) / (y[1] - 2*y[2] + y[3]) )
      table.sort( partition )
      return partition
  end,
  function(x,y,d2,d3,d4,px,py)
    local partition = {0,1}
    local p = {-d4, 3*d3, -3*d2}
    local q = {d3^2 - d2*d4, - d2*d3, d2^2}
    local inflections, doublePts = roots[2](p), roots[2](q)
    local P = multP(difP(px),difP(py))
    local monotone_partition = roots[min(#P-1,4)](P,0,1)
    for i,e in ipairs(inflections) do partition[#partition+1] = truncate_parameter(e) end
    for i,e in ipairs(doublePts) do partition[#partition+1] = truncate_parameter(e) end
    for i,e in ipairs(monotone_partition) do partition[#partition+1] = truncate_parameter(e) end
    table.sort( partition )
    return partition
  end,
  r = function(x,y,w)
    -- Find maxima
    t = {}
    t[1], t[6] = 0, 1
    t[2], t[3] = compute_rational_maxima(x[1], x[2], x[3], w)
    t[4], t[5] = compute_rational_maxima(y[1], y[2], y[3], w)
    table.sort( t )
    return t
  end
}

local function intersectLine(l1,l2) -- Finds the intersection of two lines
  local a1,b1,a2,b2 = l1[4]-l1[2], l1[1]-l1[3], l2[4]-l2[2], l2[1]-l2[3]
  local c1,c2,s1,s2 = -(a1*l1[1] + b1*l1[2]), -(a2*l2[1] + b2*l2[2]), sign(a1), sign(a2)
  local A = _M.xform(a1,b1,c1,a2,b2,c2,0,0,1)
  if not significant(a1*b2-b1*a2) then return l1[1],l1[2] end
  A = A.inverse(A)
  return A.apply(A,0,0,1)
end

-------------------------------Winding Number-----------------------------------
local winding = {}

function winding.line(c,x,y)
  if y > c.ymin and y <= c.ymax and c.l(x,y) < 0 then
    return c.s
  else return 0 end
end

function winding.quadratic(c,x,y)
  if y > c.ymin and y <= c.ymax then
    if c.l(x,y) <= 0 then
      if c.sl > 0 then return c.s
      elseif c.R(x,y)*c.sR > 0 then return c.s
      end
    else
      if c.sl > 0 and c.R(x,y)*c.sR < 0 then return c.s
      end
    end
  end
  return 0
end

function winding.cubic(c,x,y)
  if y > c.ymin and y <= c.ymax then
    if c.l(x,y) <= 0 then
      if c.sl > 0 then return c.s
      else
        if c.sT2*c.l2(x,y) > 0 and c.sT3*c.l3(x,y) > 0 then
          if c.R(x,y)*c.sR > 0 then return c.s end
        else return c.s
        end
      end
    else
      if c.sl > 0 and c.sT2*c.l2(x,y) >= 0 and c.sT3*c.l3(x,y) >= 0 then
        if c.R(x,y)*c.sR < 0 then return c.s end
      end
    end
  end
  return 0
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
  end,
  function (x_,y_,w_)
    local x, y, w = bezier_to_poly2(x_,y_,w_)
    -- This comes from the notation f10 = f1*g0 - f0*g1 used when calculating the Cailey Bézout matrix
    local f10 = {x[2]*y[1] - x[1]*y[2], w[2]*x[1] - w[1]*x[2], w[1]*y[2] - w[2]*y[1]}
    local f20 = {x[3]*y[1] - x[1]*y[3], w[3]*x[1] - w[1]*x[3], w[1]*y[3] - w[3]*y[1]}
    local f21 = {x[3]*y[2] - x[2]*y[3], w[3]*x[2] - w[2]*x[3], w[2]*y[3] - w[3]*y[2]}

    local p = prod_biP(f10,f21)
    local q = prod_biP(f20,f20)
    local R = {0,0,p[1]-q[1],1,0,p[2]-q[2],0,1,p[3]-q[3],1,1,p[4]-q[4],2,0,p[5]-q[5],0,2,p[6]-q[6]}
    return biP_to_fun(R)
  end,
  function (x_,y_)
    local x, y = bezier_to_poly3(x_,y_)
    local f10 = {x[2]*y[1] - x[1]*y[2], - x[2], y[2]}
    local f20 = {x[3]*y[1] - x[1]*y[3], - x[3], y[3]}
    local f30 = {x[4]*y[1] - x[1]*y[4], - x[4], y[4]}
    local f21, f31, f32 = x[3]*y[2] - x[2]*y[3], x[4]*y[2] - x[2]*y[4], x[4]*y[3] - x[3]*y[4]

    local p1 = prod_biP(f20,f30)
    local p2 = prod_biP(f20,f20)
    local p3 = prod_biP(f30,f30)
    f30[1] = f30[1] + f21
    local p4 = prod_biP(f10,f30)
    local p5 = prod_biP_2(f30,p3)

    local R = {0,0, f32*p4[1] + 2*f31*p1[1] - p5[1] - (f31^2)*f10[1] - f32*p2[1],
               1,0, f32*p4[2] + 2*f31*p1[2] - p5[2] - (f31^2)*f10[3] - f32*p2[2],
               0,1, f32*p4[3] + 2*f31*p1[3] - p5[3] - (f31^2)*f10[2] - f32*p2[3],
               1,1, f32*p4[4] + 2*f31*p1[4] - p5[4] - f32*p2[4],
               2,0, f32*p4[5] + 2*f31*p1[5] - p5[5] - f32*p2[5],
               0,2, f32*p4[6] + 2*f31*p1[6] - p5[6] - f32*p2[6],
               2,1, -p5[7],
               1,2, -p5[8],
               3,0, -p5[9],
               0,3, -p5[10]}
    return biP_to_fun(R)
  end
}

local implicit = {
  function (x0,y0,x1,y1)
    local l,s = implicitEq[1](x0,y0,x1,y1)
    return {x = {x0,x1}, y = {y0,y1}, ymin = min(y0,y1), ymax = max(y0,y1), l = l, s = s, type = "line"}
  end,
  function (x,y,w)
    local w = w or 1
    local l,s = implicitEq[1](x[1],y[1],x[3],y[3])
    local R = implicitEq[2](x,y,w)
    local invW = 1./w
    local x2, y2 = x[2]*invW, y[2]*invW
    local sR, sl = sign(R(x2,y2)), sign(l(x2,y2))
    return {x = x, y = y, ymin = min(y[1],y[3]), ymax = max(y[1],y[3]),
            R = R, l = l, s = s, sR = sR, sl = sl, type = "quadratic"}
  end,
  function (x,y)
    local l,s = implicitEq[1](x[1],y[1],x[4],y[4])
    local R = implicitEq[3](x,y)
    if x[1]==x[2] and x[1]==x[3] and y[1]==y[2] and y[1]==y[3] then x[2],y[2],x[3],y[3] = x[4],y[4],x[4],y[4]
    elseif x[1]==x[2] and y[1]==y[2] then x[2],y[2] = x[3],y[3]
    elseif x[4]==x[3] and y[4]==y[3] then x[3],y[3] = x[2],y[2] end

    local l2 = implicitEq[1](x[1],y[1],x[2],y[2])
    local l3 = implicitEq[1](x[3],y[3],x[4],y[4])
    local xT,yT = intersectLine({x[1],y[1],x[2],y[2]},{x[3],y[3],x[4],y[4]})
    local sR,sl,sT2,sT3 = sign(R(xT,yT)), sign(l(xT,yT)), sign(l2(x[4],y[4])), sign(l3(x[1],y[1]))

    return {x = x, y = y, ymin = min(y[1],y[4]), ymax = max(y[1],y[4]), R = R, l = l, l2 = l2, l3 = l3,
            s = s, sT2 = sT2, sT3 = sT3, sR = sR, sl = sl, type = "cubic"}
  end
}


-- prepare grid
-------------------------------------------------------------------------------------
	-------------------------------- GRID FUNCTIONS -------------------------------------
	-------------------------------------------------------------------------------------l]

	-- Contains (1) the cell coordinates and (2) the initial winding number increment
	local cell_size = 1

	--DONE AND WORKING
	local function computeGridDimension(scene) -- This can be improved to choose optimal dimensions
	    return floor(scene.height/cell_size), floor(scene.width/cell_size)
	end

	--DONE AND WORKING
	local function findCellCoord(x, y, grid) -- 0 and >n or >m will mean out of the grid
    return floor(grid.n*y*(1/grid.height) + 1), floor(grid.m*x*(1/grid.width) + 1)
	end

	-----------------------------------------------------------------------------------
  -- local function insideGrid(i,j,m,n) -- Will try to use this later
  --   return  (i > 0 and j > 0 and i < n+1 and j < m+1)
  -- end

  local function alocate_and_insert(grid,i,j,i_path,segment) -- We want to alocate space in the grid only when it is really necessary
    grid[i] = grid[i] or {}
    if grid[i][j] then
      if grid[i][j].order[#grid[i][j].order] ~= i_path then table.insert(grid[i][j].order, i_path) end
    else grid[i][j] = {order = {i_path}} end

    grid[i][j][i_path] = grid[i][j][i_path] or {initialWindingNumber = 0, segments = {}}
    table.insert(grid[i][j][i_path].segments, segment)
  end

  local function alocate_and_add_winding(grid,i,j,i_path,s)
    grid[i] = grid[i] or {}
    if grid[i][j] then
      if grid[i][j].order[#grid[i][j].order] ~= i_path then table.insert(grid[i][j].order, i_path) end
    else grid[i][j] = {order = {i_path}} end

    grid[i][j][i_path] = grid[i][j][i_path] or {initialWindingNumber = 0, segments = {}}
    grid[i][j][i_path].initialWindingNumber = grid[i][j][i_path].initialWindingNumber + s
  end

	--DONE AND WORKING
	local function walkInPath(grid, segments, i_path)
    local event_list = {}
    local n, m = grid.n, grid.m
    local x, y = {}, {}
    local xmax, xmin, ymax, ymin = 0, 0, 0, 0
    local i, j = 1, 1

    local goup = function()
      i = i+1
      event_list[#event_list+1] = {1,i,j}
    end
    local godown = function()
      event_list[#event_list+1] = {-1,i,j}
      i = i-1
    end
    local goright = function()
      j = j+1
      if y[#y] < ymax then
        alocate_and_insert(grid,i,j-1,i_path,implicit[1](x[#x],y[#y],x[#x],ymax)) --insert line going up
      end
    end
    local goleft = function()
      j = j-1
      if y[1] < ymax  then
        alocate_and_insert(grid,i,j,i_path,implicit[1](x[1],ymax,x[1],y[1])) --insert line going down
      end
    end

	  for l, segment in ipairs(segments) do
	    x, y = segment.x, segment.y
      i, j = findCellCoord(x[1], y[1], grid)

      local finali, finalj      = findCellCoord(x[#x], y[#y], grid)
      local segment_going_up    = (finali >= i)
      local segment_going_right = (finalj >= j)

      local it, Nmax = 1, 500
      while it < Nmax  do
        it = it+1

        alocate_and_insert(grid, i, j, i_path, segment)
        if (i == finali and j == finalj) then
          break --leaves the while loop and goes to another segment
        end

        xmin, xmax = (j-1)*grid.cell_width , j*grid.cell_width
  	    ymin, ymax = (i-1)*grid.cell_height, i*grid.cell_height
        local w_up_left, w_up_right, w_down_left, w_down_right = 0, 0, 0, 0

        if segment_going_up then
          w_up_left  = winding[segment.type](segment,xmin,ymax)
          w_up_right = winding[segment.type](segment,xmax,ymax)

          if (x[1] == xmin and x[#x] == xmin and y[#y] >= ymax) then
            goup()
          elseif (w_up_left ~= 0 and w_up_right == 0) or
                 (w_up_right ~= 0 and w_up_left == 0) or
                 (y[#y] == ymax and x[#x] >= xmin and x[#x] <= xmax) then
            goup()
          elseif segment_going_right and x[#x]then
            goright()
          elseif x[#x] < xmin then
            goleft()
          end
        else
          w_down_left =  winding[segment.type](segment,xmin,ymin)
          w_down_right = winding[segment.type](segment,xmax,ymin)

          if ymin == y[#y] and x[#x] >= xmin and x[#x] < xmax then
            event_list[#event_list+1] = {-1,i,j}
            break
          end
          if (y[1] == ymin and y[#y] < ymin) then
            godown()
          elseif (x[1] == xmin and x[#x] == xmin and y[#y] < ymin) then
            godown()
          elseif (w_down_left ~= 0 and w_down_right == 0) or
                 (w_down_right ~= 0 and w_down_left == 0) then
            godown()
          elseif segment_going_right then
            goright()
          elseif x[#x] < xmin then                -- if go left
            goleft()
          elseif x[#x] <= xmin then
            break
          end
        end

      end
    end
    return event_list
	end

	local function CountingSort(array, key)
		local count = {}
    local ind = {}
    if #array == 0 then return {} end
	  for j = 1, #array do
	  	i = array[j][key]
      ind[#ind+1] = i
	  	if not count[i] then count[i] = 0 end
	    count[i] = count[i] + 1
	  end
	  local total = 1
    local maxim = max(unpack(ind))
		for i = 1,maxim do
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

  local function running_sum(grid, list, i_path)
    local total = 0
    local wN,i,j = unpack(list[#list])
    for k = #list-1,1,-1 do
      local next_i, next_j = list[k][2], 1
      local same_line = (next_i == i)
      if same_line then
        next_j = list[k][3]
      end
      if wN ~= 0 then
        for l=next_j,j-1 do
          alocate_and_add_winding(grid,i,l,i_path,wN)
        end
      end
      if same_line then
        wN = wN + list[k][1]
      else
        wN = list[k][1]
      end
      i, j = list[k][2], list[k][3]
    end
    if wN ~= 0 then
      for l = 1,j-1 do
        alocate_and_add_winding(grid,i,l,i_path,wN)
      end
    end
  end

	function prepare.grid(scene)
    scene.grid = {}
	  scene.grid.n, scene.grid.m = computeGridDimension(scene) --n x m grid
    scene.grid.cell_width, scene.grid.cell_height = scene.width/scene.grid.m, scene.height/scene.grid.n
    scene.grid.width, scene.grid.height = scene.width, scene.height --This could've been passed directly

		for i, e in ipairs(scene.elements) do
			local event_list = walkInPath(scene.grid, e.shape.segments, i)
      if #event_list > 0 then
        local y_order = CountingSort(event_list, 3) -- (sort by y line)
        local x_order = CountingSort(y_order, 2)
        running_sum(scene.grid,x_order,i)
      end
		end
	end

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
  if shape.degenerated then
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
  local x = {  data[pos], data[pos+2], data[pos+4]}
  local y = {data[pos+1], data[pos+3], data[pos+5]}

  if colinear(x,y) then --If it degenerates to a line
    shape.segments[#shape.segments+1] = implicit[1](x[1],y[1],x[3],y[3])
    return shape
  end

  local S = partition[2](x,y)
  for i = 2, 4 do
      if S[i-1] ~= S[i] then
        local x_, y_ = bezier_cut[2](S[i-1], S[i], x, y)
        shape.segments[#shape.segments+1] = implicit[2](x_,y_)
      end
  end
end

function execute.cubic_segment(shape)
  local data, pos = shape.data, shape.pos
  dataTransform(shape,data,pos+2,3)

  local x = {  data[pos], data[pos+2], data[pos+4], data[pos+6]}
  local y = {data[pos+1], data[pos+3], data[pos+5], data[pos+7]}

  --Finding inflections and double points
  local px, py = bezier_to_poly3(x,y)
  local d2,d3,d4 = cross_with_e1(px,py)

  local degree = 3
  if d2==0 and d3==0 then
    if d4==0 then--Case it degenerates to a straight line
      shape.segments[#shape.segments+1] = implicit[1](x[1],y[1],x[4],y[4])
      return shape
    else --Case it degenerates to a quadratic
      degree = 2
      local xaux = (3/2)*(x[2] - x[1]) + x[1]
			local yaux = (3/2)*(y[2] - y[1]) + y[1]
      x[2],y[2],x[3],y[3] = xaux,yaux,x[4],y[4]
    end
  end

  local S = partition[degree](x,y,d2,d3,d4,px,py)
  for i = 2, #S do
      if S[i-1] ~= S[i] then
        local x_,y_ = bezier_cut[degree](S[i-1], S[i], x, y)
        shape.segments[#shape.segments+1] = implicit[degree](x_,y_)
      end
  end
end

function execute.rational_quadratic_segment(shape)
  local data, pos = shape.data, shape.pos
  data[pos+2], data[pos+3], data[pos+4] = shape.xf : apply(data[pos+2], data[pos+3], data[pos+4])
  data[pos+5], data[pos+6] = shape.xf : apply(data[pos+5], data[pos+6])

  local x = {  data[pos], data[pos+2], data[pos+5]}
  local y = {data[pos+1], data[pos+3], data[pos+6]}
  local w = data[pos+4]
  local invW = 1./w
  local x_p = {x[1], x[2]*invW, x[3]}
  local y_p = {y[1], y[2]*invW, y[3]}
  shape.xmin = min(shape.xmin, unpack(x_p))
  shape.ymin = min(shape.ymin, unpack(y_p))
  shape.xmax = max(shape.xmax, unpack(x_p))
  shape.ymax = max(shape.ymax, unpack(y_p))

  if colinear(x_p,y_p) then --If it degenerates to a line
    shape.segments[#shape.segments+1] = implicit[1](x[1],y[1],x[3],y[3])
    return shape
  end

  local S = partition["r"](x,y,w)
  for i = 2, 6 do
      if S[i-1] ~= S[i] then
          local x_, y_, w_ = bezier_cut["r"](S[i-1], S[i], x, y, w)
          shape.segments[#shape.segments+1] = implicit[2](x_,y_,w_)
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
  return shape
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
  return shape
end

function prepare.circle(shape)
  -- it is formed by 4 arcs covering each quarter of the unit circle (Using Four Vertex Theorem[Manfredo C.] )
  local s = sqrt(2)/2          -- sin(pi/2)
  local cx, cy, r, xf = shape.cx, shape.cy, shape.r, shape.xf
  shape = _M.path{
      _M.M,  1,  0,
      _M.R,  s,  s,  s, 0, 1,
      _M.R, -s,  s,  s,-1, 0,
      _M.R, -s, -s,  s, 0,-1,
      _M.R,  s, -s,  s, 1, 0,		-- colocar em baixo os parametros do circulo
      _M.Z}:scale(r, r):translate(cx, cy)
  shape.xf = xf*shape.xf
  return prepare.path(shape)
end
--------------------------------------------------------------------------------
--------------------------Paint auxiliares--------------------------------------
local spread = {
  pad = function (t)
    return min(1,max(0,t))
  end,
  ["repeat"] = function (t)
    return t-floor(t)
  end,
  reflect = function (t)
    return 2*math.abs(t*.5 - floor(t*.5 + .5))
  end,
  transparent = function (t) return t end
}

local function interpolate(t1,t2,t,c1,c2)
  return mult(1/(t2-t1),addM(mult(t2-t,c1),mult(t-t1,c2)))
end

local function ramp(r,t)
  for i=1,(#r-3),2 do
    if t >= r[i] and t < r[i+2] then
      return interpolate(r[i],r[i+2],t,r[i+1],r[i+3])
    end
  end
  if r.spread == "transparent" then return {0,0,0,0}
  elseif t < r[1] then return r[2]
  elseif t >= r[#r-1] then return r[#r]
  end
  return {0,0,0,0}
end

local function ramp_texture(img,p)
  local w, h = img.width, img.height
  local x, y = p[1]*w + 1/2, p[2]*h + 1/2
  local x1,y1,x2,y2,x3,y3,x4,y4 = floor(x), ceil(y), ceil(x), ceil(y), floor(x), floor(y), ceil(x), floor(y)
  local a1,a2,a3,a4 = (x4-x)*(y-y4), (x-x3)*(y-y3), (x2-x)*(y2-y), (x-x1)*(y1-y)
  local c1,c2,c3,c4 = {img:get(max(x1,1),min(y1,h))}, {img:get(min(x2,w),min(y2,h))}, {img:get(max(x3,1),max(y3,1))}, {img:get(min(x4,w),max(y4,1))}
  return addM(addM(mult(a1,c1),mult(a2,c2)),addM(mult(a3,c3),mult(a4,c4)))
end



--------------------------Preparing Paint---------------------------------------
function prepare.solid(paint)
  paint.data[4] = paint.data[4]*paint.opacity
  paint.color = function (self,x,y)
    return self.data
  end
end

function prepare.lineargradient(paint)
  local data = paint.data
  for i=1,#data.ramp-1,2 do
    data.ramp[i+1][4] = data.ramp[i+1][4]*paint.opacity
  end
  local wrapping = spread[data.ramp.spread]
  local function lmap(p)
    local v = data.p2 - data.p1
    return _M.dot(p-data.p1,v)/(_M.dot(v,v)) -- a lot of divisions :/
  end
  paint.color = function (self,x,y)
    return ramp(data.ramp,wrapping(lmap(self.xf*_M.vector(x,y))))
  end
end

function prepare.radialgradient(paint)
  local data = paint.data
  for i=1,#data.ramp-1,2 do
    data.ramp[i+1][4] = data.ramp[i+1][4]*paint.opacity
  end
  local r,c,f = data.radius, data.center, data.focus
  assert(r > 0, "wrong parameters")
  local wrapping = spread[data.ramp.spread]
  local v = f-c
  local d_to_v = _M.dot(v,v)
  if d_to_v > r^2 then
    r = c + (r/math.sqrt(d_to_v))*v
  end
  local function rmap(p)
    local u = p-f
    local d1 = _M.dot(u,u)
    if d1 == 0 then return 0 end
    local p = {d_to_v-r^2,2*_M.dot(u,v),d1}
    local t = roots[2](p)
    if t[1] > 0 then return 1/t[1]
    elseif t[2] > 0 then return 1/t[2]
    else return 1
    end
  end
  paint.color = function (self,x,y)
    return ramp(data.ramp,wrapping(rmap(self.xf*_M.vector(x,y))))
  end
end

function prepare.texture(paint)
  local img = paint.data.image
  local wrapping = function(x,y) return _M.vector(spread[paint.data.spread](x), spread[paint.data.spread](y)) end
  paint.color = function (self,x,y)
    return ramp_texture(img,wrapping(unpack(self.xf*_M.vector(x,y))))
  end
end

--------------------------------------------------------------------------------
-- prepare scene for sampling and return modified scene
local function preparescene(scene)
  for i, element in ipairs(scene.elements) do
    element.shape.xf = scene.xf * element.shape.xf
    element.shape = prepare[element.shape.type](element.shape)
    element.paint.xf = element.paint.xf:inverse() * scene.xf:inverse()
    prepare[element.paint.type](element.paint)
  end
  prepare.grid(scene)
  return scene
end

local gamma_factor = 2.2
local n_samples = 16
local sampling_pattern = blue[n_samples]

local function gamma(c)
    local f = 1/gamma_factor
    return {c[1]^f, c[2]^f, c[3]^f, c[4]^f}
end

local function ungamma(c)
    return {c[1]^gamma_factor,c[2]^gamma_factor,c[3]^gamma_factor,c[4]^gamma_factor}
end

local function sample_point(scene, x, y)
  local color = {1,1,1,1}
  local k,l = findCellCoord(x,y,scene.grid)
  if scene.grid[k] == nil then return color end
  if scene.grid[k][l] then
    for ind,i in ipairs(scene.grid[k][l].order) do
      local cell = scene.grid[k][l][i]
      local shape, paint = scene.elements[i].shape, scene.elements[i].paint
      if inside(shape,x,y) then
        local wN = cell.initialWindingNumber + windingNumber(cell.segments,x,y)
        if scene.elements[i].type == "fill" and wN ~= 0 then
          color = blend(paint:color(x,y),color)
        elseif scene.elements[i].type == "eofill" and wN % 2 == 1 then
          color = blend(paint:color(x,y),color)
        end
      end
    end
  end
  return color
end

-- sample scene at x,y and return r,g,b,a
local function sample(scene, x, y)
  local sum = {0,0,0,0}
  for i = 1,#sampling_pattern,2 do
    local dx, dy = sampling_pattern[i], sampling_pattern[i+1]
    sum = addM(sum, ungamma(sample_point(scene, x + dx, y + dy)))
  end
  return unpack(gamma(mult(1/n_samples,sum)))
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

function _M.render(scene, viewport, file, arguments)
local time = chronos.chronos()
local preprocess = chronos.chronos()
local pattern = blue[1]
local maxdepth = MAX_DEPTH
local minsegs = MIN_SEGS
local scenetree = false
local prefix = "none"
if #arguments > 0 then stderr("driver arguments:\n") end
for i, argument in ipairs(arguments) do
    stderr("  %d: %s\n", i, argument)
end
local tx, ty = 0., 0.
local p = nil
-- list of supported options
local options = {
    { "^%-prefix:(.*)$", function(n)
        if not n then return false end
        prefix = n
        return true
    end },
    { "^(%-maxdepth:(%d+)(.*))$", function(all, n, e)
        if not n then return false end
        assert(e == "", "invalid option " .. all)
        maxdepth = assert(tonumber(n), "invalid option " .. all)
        assert(maxdepth >= 1, "invalid option " .. all)
        return true
    end },
    { "^(%-tx:(%-?%d+)(.*))$", function(all, n, e)
        if not n then return false end
        assert(e == "", "trail invalid option " .. all)
        tx = assert(tonumber(n), "number invalid option " .. all)
        return true
    end },
    { "^(%-ty:(%-?%d+)(.*))$", function(all, n, e)
        if not n then return false end
        assert(e == "", "trail invalid option " .. all)
        ty = assert(tonumber(n), "number invalid option " .. all)
        return true
    end },
    { "^(%-p:(%d+)(.*))$", function(all, n, e)
        if not n then return false end
        assert(e == "", "trail invalid option " .. all)
        p = assert(tonumber(n), "number invalid option " .. all)
        assert(p >= 1 and p <= #scene.elements, "path out of range")
        return true
    end },
    { "^(%-pattern:(%d+)(.*))$", function(all, n, e)
        if not n then return false end
        assert(e == "", "trail invalid option " .. all)
        n = assert(tonumber(n), "number invalid option " .. all)
        assert(blue[n], "non exist invalid option " .. all)
        pattern = blue[n]
        return true
    end },
    { "^(%-minsegs:(%d+)(.*))$", function(all, n, e)
        if not n then return false end
        assert(e == "", "invalid option " .. all)
        minsegs = assert(tonumber(n), "invalid option " .. all)
        return true
    end },
    { "^%-verbose$", function(all)
        if not all then return false end
        verbose = true
        return true
    end },
    { ".*", function(all)
        error("unrecognized option " .. all)
    end }
}
-- process options
for i, argument in ipairs(arguments) do
    for j, option in ipairs(options) do
        if option[2](argument:match(option[1])) then
            break
        end
    end
end
if p then scene.elements = { scene.elements[p] } end
-- potentially apply translation to scene for debugging
scene = scene:translate(tx, ty)
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
