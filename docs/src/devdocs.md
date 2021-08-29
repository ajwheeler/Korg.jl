This page contains info for people interested in contributing code to Korg.  If you 
have questions, do not hesitate to ask.

## code guidelines
- Try to be explicit about units throughout the code, particularly when not using CGS.
- Whenever possible, calculations should be precise up to a factor of $$10^{-3}$$.  When it's easy and inexpensive, they should be precise to $$10^{-5}$$ or better.  
- Ensure types are generic enough to support dual numbers and autodifferentiation. 
- Limit lines to 100 characters.

## Complete API
Here are all the documented methods in Korg.

```@autodocs
Modules = [Korg, Korg.ContinuumOpacity]
```
