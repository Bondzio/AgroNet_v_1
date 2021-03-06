USE [DEAQ 20141007]
GO
/****** Object:  StoredProcedure [dbo].[plm_spGetPProductCrops]    Script Date: 09/09/2015 18:22:26 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

/*
	Author:			Miguel Ramírez/ Nalleli López			 
	Object:			dbo.[plm_spGetPProductCrops]
	
	Company:		PLM.
	
	EXECUTE dbo.[plm_spGetPProductCrops]  @ProductId = 23333

*/


CREATE PROCEDURE [dbo].[plm_spGetPProductCrops]

(
	@ProductId int
)
AS
BEGIN

		select * from Crops a
		left join ProductCrops ps on a.CropId = ps.CropId
								and ProductId = @ProductId
		where ps.ProductId is null order by 2
END