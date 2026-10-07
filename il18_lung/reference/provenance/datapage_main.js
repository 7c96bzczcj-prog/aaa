$(document).ready(function(){



  var dataTable = $('#dataPage_table').DataTable({
        "paging":   false,
        "ordering": false,
        "info":     false,
        "searching": true,
        "bPaginate": false,
        "bInfo": false,
        "bFilter": false,
        "processing": true,
        "serverSide": false ,
        "searchHighlight": true,
        "ajax": {
          "url": 'dataPage_get_data.php',
          // "url": 'metainfo_table_get_data.php',
          "type": 'POST',
          // "data": {
          //   "genesymbol": datum
          // }
        },

        columns: [
          {
            "data": "Name"
          },
          {
            "data": "DataTypes"
          },
          {
            "data":  "GeoAccesion",
            "render":function(data,type,row){
              return  '<a href ="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+ row.GeoAccession +'" target= "_blank">'+row.GeoAccession+ '</a>';
            }

          },
          {
            "data": "Description"
          },
          {
            "data": "CellTypes"
          },

          {
            "data":  "PubMed",
            "render": function(data, type, row){
              return '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed +'" target= "_blank">'+row.PubMed+ '</a>'+ ' ' + '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed2 +'" target= "_blank">'+row.PubMed2+ '</a>' + ' '
              + '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed3 +'" target= "_blank">'+row.PubMed3+ '</a>'+ ' ' + '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed4 +'" target= "_blank">'+row.PubMed4+ '</a>' + ' ' +
              '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed5 +'" target= "_blank">'+row.PubMed5+ '</a>' + ' ' + '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed6 +'" target= "_blank">'+row.PubMed6+ '</a>' + ' ' +
              '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed7 +'" target= "_blank">'+row.PubMed7+ '</a>' + ' ' + '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed8 +'" target= "_blank">'+row.PubMed8+ '</a>' + ' ' +
              '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed9 +'" target= "_blank">'+row.PubMed9+ '</a>' + ' ' + '<a href ="https://www.ncbi.nlm.nih.gov/pubmed/'+ row.PubMed10 +'" target= "_blank">'+row.PubMed10+ '</a>'
            }
          },
          {
            "data":  function(data, type, row, meta){

              if(data.NormalizedData == "Not Applicable"){
              	  if(data.GeoAccession == 'GSE100738'){
              	  	  return '<a href="https://sharehost.hms.harvard.edu/immgen/ImmGenATAC18_AllOCRsInfo.csv" target= "_blank">Processed ATAC-seq data and called peaks</a>';
              	  }
              	  else{
              	  	  return "Not Applicable";
                  }
              }
              else{
                return '<a href="https://sharehost.hms.harvard.edu/immgen/'+ data.GeoAccession +'/'+data.NormalizedData+'.csv" target= "_blank">'+data.NormalizedName +'</a>';
              }

            }
            // "data": "NormalizingData",
            // "render":function(data, type, row){
            //   return '<a href ="https://sharehost.hms.harvard.edu/immgen/'+row.GeoAccesion+'/'+ row.NormalizingData +'.csv" target= "_blank">'+row.NormalizingData+ '</a>'
            // }
          },{
            "data": function(data, type, row, meta){

              if(data.RawGeneCount == "Not Applicable"){
                return "Not Applicable";
              }
              else{
                  return '<a href="https://sharehost.hms.harvard.edu/immgen/'+ data.GeoAccession +'/'+data.RawGeneCount+'.csv" target= "_blank">'+data.RawName +'</a>';
              }
            }
          }


        ],

        "columnDefs": [
            {
                "width": 200,
                "targets": [ 1 ],
                "visible": false,
                "searchable": true
            }
        ],

        initComplete: function () {
          this.api().columns(1).every( function () {
              var column = this;
              var select = $('<select><option value="">All</option></select>')
                  .appendTo('#drop_menu' )
                  .on( 'change', function () {
                      var val = $.fn.dataTable.util.escapeRegex(
                          $(this).val()
                      );

                      column
                          .search( val ? '^'+val+'$' : '', true, false )
                          .draw();
                  } );

              column.data().unique().sort().each( function ( d, j ) {
              		  if(d == 'z_Other'){
              		  	  select.append( '<option value="z_Other">Other</option>' )
              		  } else {
              		  	  select.append( '<option value="'+d+'">'+d+'</option>' )
                  	}
              } );
          } );
      },

        "pagingType": "simple"
      });



      dataTable.columns().every( function () {
      var column = this;

      $( 'input', this.header() ).on( 'keyup change', function () {
          column
              .search( this.value )
              .draw();


            // column.highlight(dataTable.search());
          // column.unhighlight();
          // column.highlight( column.search() );
      } );
    } );


      $("#myInputText").on("keyup", function() {
          // document.getElementById("dataPage_panel").style.display = "block";
          dataTable.search($(this).val()).draw();
          // $(this).toggle($(this).text().toLowerCase().indexOf(value) > -1)
      });


});
