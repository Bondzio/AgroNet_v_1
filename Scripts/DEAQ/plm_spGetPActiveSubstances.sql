USE [DEAQ 20141007]
GO
/****** Object:  StoredProcedure [dbo].[plm_spGetPActiveSubstances]    Script Date: 09/09/2015 18:15:48 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

/*
	Author:			Miguel Ramírez/ Nalleli López			 
	Object:			dbo.[plm_spGetPActiveSubstances]
	
	Company:		PLM.
	
	EXECUTE dbo.[plm_spGetPActiveSubstances]  @ProductId = 23333

*/


CREATE PROCEDURE [dbo].[plm_spGetPActiveSubstances]
(
	@ProductId int
)
AS
BEGIN
	select * from ActiveSubstances a
	left join ProductSubstances ps on a.ActiveSubstanceId = ps.ActiveSubstanceId
								and ProductId = @ProductId
	where ps.ProductId is null order by 2

END