--This function is to find the remainder of two integers x,y, where y is not 0.
--Input: x,y : two integers
--Output: the absoute value of remainder of x,y, whose absoute value is smaller than half of absoute value of y.
local remainder_half = function(x, y)
	if x % y == 0 then
		return 0
	elseif y % 2 ~= 0 then
		if math.abs(x % y) <= math.abs(y) / 2 then
			return math.abs(x % y)
		else
			return math.abs(x % y - y)
		end
	else
		if math.abs(x % y) == math.abs(y) / 2 then
			return math.abs(x % y)
		else
			if math.abs(x % y) < math.abs(y / 2) then
				return math.abs(x % y)
			else
				return math.abs(x % y - y)
			end
		end
	end
end
--This function is to find the maximum common factor of two integels x,y, where x is not 0 or y is not 0.
--Input: x,y: two integels
--Output: gcd(x,y): the maximum common factor of x,y
local function gcd(x, y)
	if math.abs(x) > math.abs(y) then
		local a = x
		x = y
		y = a
	end
	if x == 0 then
		return math.abs(y)
	else
		return gcd(remainder_half(y, x), x)
	end
end
