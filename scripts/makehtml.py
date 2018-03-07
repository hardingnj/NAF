import pandas as pd
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models.widgets import DataTable, TableColumn, Select
from bokeh.plotting import save, output_file
from bokeh.models.widgets import NumberFormatter

df = pd.read_csv(snakemake.input.csv, index_col=0)

source = ColumnDataSource(df)
original_source = ColumnDataSource(df)

numform = NumberFormatter(format='0.00')

strcolumns = [TableColumn(field=x, title=x) for x in ["coach", "race", "nation"]]
numcolumns = [TableColumn(field=x,
                          title=x, 
                          formatter=numform) for x in ["mu", "phi", "value"]]

data_table = DataTable(source=source, columns=strcolumns + numcolumns)

data_table.height = 1000
data_table.width = 500

# callback code to be used by all the filter widgets
# requires (source, original_source, race_select_obj,target_object)
combined_callback_code = """
var data = source.get('data');
var original_data = original_source.get('data');
var race = race_select_obj.get('value');
console.log("race: " + race);

for (var key in original_data) {
    data[key] = [];
    for (var i = 0; i < original_data['race'].length; ++i) {
        if (race === "ALL" || original_data['race'][i] === race) {
            data[key].push(original_data[key][i]);
        }
    }
}

target_obj.trigger('change');
source.trigger('change');

"""

# define the filter widgets, without callbacks for now
race_list = ['ALL'] + df['race'].unique().tolist()
race_select = Select(title="Race:", value=race_list[0], options=race_list)

# year_list = ['ALL'] + df['year'].unique().tolist()
# year_select = Select(title="Year:", value=year_list[0], options=year_list)
# now define the callback objects now that the filter widgets exist
generic_callback = CustomJS(
    args=dict(source=source, 
              original_source=original_source, 
              race_select_obj=race_select, 
              target_obj=data_table),
    code=combined_callback_code
)

# finally, connect the callbacks to the filter widgets
race_select.js_on_change('value', generic_callback)

#year_select.js_on_change('value', generic_callback)
p = column(race_select, data_table)

output_file('output/rankings.html')

save(p)
